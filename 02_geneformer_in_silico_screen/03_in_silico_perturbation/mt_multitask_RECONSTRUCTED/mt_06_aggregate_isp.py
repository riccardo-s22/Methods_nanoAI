"""Aggregate the MT in-silico-perturbation output into per-gene statistics.

Runs Geneformer's InSilicoPerturberStats for the significance test, then pools
the per-nucleus shifts into per-gene statistics, applies the 20%-detection
eligibility rule, and - this is what distinguishes the MT aggregation from the
FT one - writes a separate table per annotated cell type as well as the pooled
table.

Methods ("In silico perturbation statistics"):
    InSilicoPerturberStats(mode="goal_state_shift", genes_perturbed="all",
    combos=0). No external null was supplied, so the random-gene background is
    built internally: shifts of all perturbed genes across all nuclei are pooled
    into one distribution, down-sampled to 100,000 values with a fixed seed, and
    each gene's per-nucleus shifts are compared against it by two-sided Wilcoxon
    rank-sum test with Benjamini-Hochberg FDR correction across all tested genes.
    Effect size = mean shift toward the goal state (Shift_to_goal_end).

    Eligibility: genes detected (nonzero counts) in >= 20% of nuclei of either
    the TBI or the control group. For the MT model, pooled tables use
    population-level detection and per-cell-type tables use within-cell-type
    detection, so each table resets the bar against its own population.

    "Shifts were stratified by annotated cell type, and per-gene statistics were
    computed separately for each cell type with the same test, safeguards and
    thresholds as for the FT model; excitatory neurons are reported in the main
    text."

DONOR-CONSISTENCY SAFEGUARD IS UNAVAILABLE FOR THIS COHORT. The split holds out
one donor per condition, so the 1,846 held-out TBI nuclei all come from a single
donor: `n_donors_contributing` is 1 and `donor_consistency_score` is NaN for
every gene. Do not filter on or report it for the MT screen. The remaining
safeguards (FDR via Wilcoxon vs random, n_cells_perturbed coverage,
percent-of-start-cells coverage, mean-shift magnitude) are unaffected.

=============================================================================
Rewritten from the Methods after the original MT code was lost. It documents the
procedure but was not used to generate the published results.
=============================================================================
The statistics (including the 20% eligibility rule) come from the original
`perturb_safeguards.py` in the FT folder, which is imported here; only the driver
around it was rewritten.

UNRECORDED: the exact raw-pickle filename pattern the MT runs produced. The
regex in `mt_parse_run()` follows the FT convention; adjust it to match the
actual output of mt_05_in_silico_perturbation.py.

Original file name: aggregate_garza_isp.py
=============================================================================
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

# perturb_safeguards.py and biological_gene_classes.py survive in the FT folder.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "ft_excitatory"))

import perturb_safeguards as safeguards  # noqa: E402

DETECTION_THRESHOLD = 0.20
STATE_KEY = "disease_binary"
CELL_TYPE_KEY = "cell_type"

ARM_STATES = {
    "rescue_delete": ("TBI", "control"),
    "rescue_overexpress": ("TBI", "control"),
    "induction_delete": ("control", "TBI"),
    "induction_overexpress": ("control", "TBI"),
}

# UNRECORDED - see header. Mirrors the FT naming convention.
MT_RUN_RE = re.compile(
    r"^in_silico_(?:delete|overexpress)_garza_mt_\w+?_chunk(\d+)"
    r"_dict_cell_embs_(\d+)batch(-?\d+)_raw\.pickle$"
)


def mt_parse_run(name: str):
    match = MT_RUN_RE.match(name)
    if match is None:
        return None
    chunk_id, cell_in_chunk, batch = match.groups()
    return int(chunk_id), int(cell_in_chunk), int(batch)


def run_isp_stats(chunks_dir: Path, out_dir: Path, start_state: str, goal_state: str) -> Path:
    """Geneformer's own significance test; writes one CSV into out_dir."""
    from geneformer import InSilicoPerturberStats

    stats = InSilicoPerturberStats(
        mode="goal_state_shift",
        genes_perturbed="all",
        combos=0,
        cell_states_to_model={
            "state_key": STATE_KEY,
            "start_state": start_state,
            "goal_state": goal_state,
            "alt_states": [],
        },
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    stats.get_stats(str(chunks_dir), None, str(out_dir), "isp_stats")
    produced = sorted(out_dir.glob("isp_stats*.csv"))
    if not produced:
        raise FileNotFoundError(f"InSilicoPerturberStats wrote no CSV into {out_dir}")
    return produced[0]


def aggregate_one_population(
    shifts: list[tuple[int, int, float]],
    donor_per_cell: list[str],
    n_start_cells: int,
    token_to_ensembl: dict[int, str],
    token_to_symbol: dict[int, str],
    isp_csv: Path,
    det_start: dict,
    det_goal: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-gene table for one population (pooled, or one cell type)."""
    per_gene = safeguards.aggregate_per_gene(
        iter(shifts),
        donor_per_cell=donor_per_cell,
        total_tbi_cells=n_start_cells,
        token_to_ensembl=token_to_ensembl,
        token_to_symbol=token_to_symbol,
    )
    per_gene = safeguards.merge_with_isp_stats(per_gene, isp_csv)
    annotated, filtered = safeguards.apply_eligibility_filter(
        per_gene, det_start=det_start, det_goal=det_goal, thresh=DETECTION_THRESHOLD
    )
    return safeguards.rank_and_order(annotated), safeguards.rank_and_order(filtered)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", required=True, choices=sorted(ARM_STATES))
    ap.add_argument("--chunks-dir", required=True, help="raw ISP pickles for this arm")
    ap.add_argument(
        "--tokenized-heldout",
        required=True,
        help="tokenized held-out-donor .dataset that was perturbed",
    )
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--cells-per-chunk", type=int, default=56)
    ap.add_argument("--skip-isp-stats", action="store_true")
    args = ap.parse_args()

    from datasets import load_from_disk

    start_state, goal_state = ARM_STATES[args.arm]
    chunks_dir, out_dir = Path(args.chunks_dir), Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Geneformer significance test -------------------------------------
    isp_csv = (
        sorted(out_dir.glob("isp_stats*.csv"))[0]
        if args.skip_isp_stats
        else run_isp_stats(chunks_dir, out_dir, start_state, goal_state)
    )
    print(f"InSilicoPerturberStats table: {isp_csv}")

    # --- the perturbed population ----------------------------------------
    ds = load_from_disk(args.tokenized_heldout)
    is_start = [s == start_state for s in ds[STATE_KEY]]
    start_rows = [i for i, keep in enumerate(is_start) if keep]
    cell_types = [ds[CELL_TYPE_KEY][i] for i in start_rows]
    donors = [ds["donor"][i] for i in start_rows]
    n_start_cells = len(start_rows)
    print(f"{args.arm}: {n_start_cells:,} start-state ({start_state}) nuclei perturbed")
    print(f"  donors contributing: {sorted(set(donors))}")

    token_to_ensembl, token_to_symbol = safeguards.load_token_maps()

    # --- pool the per-nucleus shifts --------------------------------------
    spec = safeguards.RunSpec(
        chunks_dir=chunks_dir,
        parse_fn=mt_parse_run,
        cells_per_chunk=args.cells_per_chunk,
        total_tbi_cells=n_start_cells,
        goal_state=goal_state,
    )
    shifts = list(safeguards.iter_per_cell_shifts(spec))
    print(f"  pooled {len(shifts):,} (nucleus, gene) shift entries")

    # --- detection fractions ---------------------------------------------
    det_pooled, det_byct = safeguards.compute_detection_fractions(
        args.tokenized_heldout,
        state_key=STATE_KEY,
        states=[start_state, goal_state],
        cell_type_key=CELL_TYPE_KEY,
    )

    # --- pooled table -----------------------------------------------------
    all_genes, eligible = aggregate_one_population(
        shifts,
        donor_per_cell=donors,
        n_start_cells=n_start_cells,
        token_to_ensembl=token_to_ensembl,
        token_to_symbol=token_to_symbol,
        isp_csv=isp_csv,
        det_start=det_pooled[start_state],
        det_goal=det_pooled[goal_state],
    )
    all_genes.to_csv(out_dir / f"{args.arm}_unfiltered_all_gene_ranked.tsv", sep="\t", index=False)
    eligible.to_csv(
        out_dir / f"{args.arm}_filtered_high_confidence_gene_ranked.tsv", sep="\t", index=False
    )
    print(f"  pooled: {len(all_genes):,} genes tested, {len(eligible):,} eligible")

    # --- per-cell-type tables (the MT-specific step) ----------------------
    by_ct_dir = out_dir / "by_cell_type"
    by_ct_dir.mkdir(exist_ok=True)
    for cell_type in sorted(set(cell_types)):
        rows_in_ct = {i for i, ct in enumerate(cell_types) if ct == cell_type}
        ct_shifts = [s for s in shifts if s[0] in rows_in_ct]
        if not ct_shifts:
            continue
        ct_all, ct_eligible = aggregate_one_population(
            ct_shifts,
            donor_per_cell=donors,
            n_start_cells=len(rows_in_ct),
            token_to_ensembl=token_to_ensembl,
            token_to_symbol=token_to_symbol,
            isp_csv=isp_csv,
            det_start=det_byct[(cell_type, start_state)],
            det_goal=det_byct[(cell_type, goal_state)],
        )
        ct_eligible.to_csv(
            by_ct_dir / f"{args.arm}_{cell_type}_filtered_high_confidence_gene_ranked.tsv",
            sep="\t",
            index=False,
        )
        print(
            f"  {cell_type:>22s}: {len(rows_in_ct):>5,} nuclei, "
            f"{len(ct_all):>6,} tested, {len(ct_eligible):>6,} eligible"
        )

    print(f"\nwrote tables to {out_dir}")


if __name__ == "__main__":
    main()
