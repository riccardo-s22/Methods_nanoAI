#!/usr/bin/env python
"""
Aggregate fine-tuned control→TBI InSilicoPerturber outputs into per-gene
cosine-shift statistics with coverage safeguards.

Mirrors `05_get_perturb_stats_paperQC.py` but:
  • starts from CONTROL cells — `start_state="control"`, `goal_state="TBI"`
  • supports both `delete` and `overexpress` perturbation modes via --mode
  • supports multiple cell types via --celltype (see finetuned_celltype_config.py)
  • reads / writes mode-specific subdirs under each cell type's run_root:
        <run_root>/chunks_<mode>/
        <run_root>/stats_<mode>/

Statistical core is unchanged: Geneformer `InSilicoPerturberStats(mode="goal_state_shift",
genes_perturbed="all", null_dist_data_directory=None)` — Wilcoxon rank-sum of each gene's
per-cell shifts vs. the pooled all-genes random distribution (Theodoris 2023 Fig 6a).

Usage:
  python scripts/05_get_perturb_stats_finetuned.py --mode delete      [--celltype excit|oligo] [--allow-partial] [--skip-isp]
  python scripts/05_get_perturb_stats_finetuned.py --mode overexpress [--celltype excit|oligo] [--allow-partial] [--skip-isp]
"""
import argparse
import sys
from pathlib import Path

# Local helpers
sys.path.insert(0, str(Path(__file__).resolve().parent))
# finetuned_celltype_config.py sits in 02_model_finetuning/ft_excitatory in this
# repository (the scripts were written for one flat scripts/ folder).
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "02_model_finetuning" / "ft_excitatory"))
from finetuned_celltype_config import CELLTYPE_CONFIG, get_config
from perturb_safeguards import (
    RunSpec,
    aggregate_per_gene,
    apply_high_confidence_filter,
    iter_per_cell_shifts,
    load_donor_map,
    load_token_maps,
    make_parse_full_run,
    merge_with_isp_stats,
    rank_and_order,
)

# Same coverage thresholds as the zero-shot pipeline
MIN_CELLS_PERTURBED = 25
MIN_PCT_CELLS       = 2.0

parser = argparse.ArgumentParser()
parser.add_argument("--mode", choices=["delete", "overexpress"], required=True,
                    help="Perturbation type to aggregate")
parser.add_argument("--celltype", choices=sorted(CELLTYPE_CONFIG.keys()),
                    default="excit",
                    help="Cell type config from finetuned_celltype_config.py (default: excit)")
parser.add_argument("--allow-partial", action="store_true",
                    help="Run even if fewer than expected chunks are present")
parser.add_argument("--skip-isp", action="store_true",
                    help="Skip re-running InSilicoPerturberStats; expect existing CSV")
args = parser.parse_args()

cfg = get_config(args.celltype)
RUN_ROOT            = cfg["run_root"]
TOKENIZED           = cfg["tokenized"]
CELLTYPE_PREFIX     = cfg["celltype_prefix"]
EXPECTED_CHUNKS     = cfg["expected_chunks"]
CELLS_PER_CHUNK     = cfg["cells_per_chunk"]
TOTAL_START_CELLS   = cfg["total_start_cells"]
START_STATE         = cfg["start_state"]
GOAL_STATE          = cfg["goal_state"]

# Direction-aware cell states (perturbing START → toward GOAL centroid).
CELL_STATES = {
    "state_key":   "condition",
    "start_state": START_STATE,
    "goal_state":  GOAL_STATE,
    "alt_states":  [],
}

CHUNKS_DIR = RUN_ROOT / f"chunks_{args.mode}"
STATS_DIR  = RUN_ROOT / f"stats_{args.mode}"
OUTPUT_PREFIX = f"{CELLTYPE_PREFIX}_{args.mode}"

print(f"[config] celltype={args.celltype}  mode={args.mode}")
print(f"[config] direction: {START_STATE} → {GOAL_STATE}")
print(f"[config] run_root={RUN_ROOT}")
print(f"[config] tokenized={TOKENIZED.name}")
print(f"[config] prefix={CELLTYPE_PREFIX}  expected_chunks={EXPECTED_CHUNKS}  cells_per_chunk={CELLS_PER_CHUNK}  total_start_cells={TOTAL_START_CELLS}")

parse_chunks = make_parse_full_run(CELLTYPE_PREFIX, isp_mode=args.mode)

# ── preflight ──────────────────────────────────────────────────────────────
pickles = sorted(CHUNKS_DIR.glob("*_raw.pickle"))
print(f"[{args.mode}] Found {len(pickles):,} pickle file(s) in {CHUNKS_DIR}")

if len(pickles) == 0:
    sys.exit(
        f"ERROR: No _raw.pickle files found in {CHUNKS_DIR}.\n"
        f"Ensure the Slurm array job (run_perturb_array_finetuned_control_to_tbi_{args.mode}.sh) has completed."
    )

cells_seen: set[int] = set()
for p in pickles:
    parsed = parse_chunks(p.name)
    if parsed is None:
        continue
    chunk_id, cell_in_chunk, _ = parsed
    cells_seen.add(chunk_id * CELLS_PER_CHUNK + cell_in_chunk)
print(f"Unique {START_STATE} cells with pickles: {len(cells_seen):,} / {TOTAL_START_CELLS:,} expected")

if len(cells_seen) < TOTAL_START_CELLS:
    msg = (
        f"WARNING: {len(cells_seen)}/{TOTAL_START_CELLS} {START_STATE} cells have pickles. "
        "Statistics will cover a subset."
    )
    if args.allow_partial:
        print(msg)
    else:
        sys.exit(f"{msg}\nRun with --allow-partial to aggregate anyway.")

STATS_DIR.mkdir(parents=True, exist_ok=True)
raw_csv = STATS_DIR / f"{OUTPUT_PREFIX}.csv"

# ── 1. run InSilicoPerturberStats (FDR / pval / Sig) ───────────────────────
if args.skip_isp and raw_csv.exists():
    print(f"\n[skip-isp] Using existing InSilicoPerturberStats CSV: {raw_csv}")
else:
    from geneformer import InSilicoPerturberStats

    print(f"\nRunning InSilicoPerturberStats (mode=goal_state_shift) ...")
    ispstats = InSilicoPerturberStats(
        mode="goal_state_shift",
        genes_perturbed="all",
        combos=0,
        anchor_gene=None,
        cell_states_to_model=CELL_STATES,
    )
    ispstats.get_stats(
        input_data_directory=str(CHUNKS_DIR),
        null_dist_data_directory=None,
        output_directory=str(STATS_DIR),
        output_prefix=OUTPUT_PREFIX,
    )
    if not raw_csv.exists():
        sys.exit(f"ERROR: Expected output not found at {raw_csv}")
    print(f"  Wrote: {raw_csv}")

# ── 2. compute per-cell safeguards from raw pickles ────────────────────────
print(f"\nLoading {START_STATE} donor mapping from {TOKENIZED.name} ...")
donor_per_cell, _sample_per_cell = load_donor_map(TOKENIZED, condition=START_STATE)
print(f"  {START_STATE} cells in dataset: {len(donor_per_cell):,}  unique donors: {len(set(donor_per_cell))}")

print("Loading Geneformer token → ensembl / symbol maps ...")
token_to_ensembl, token_to_symbol = load_token_maps()

spec = RunSpec(
    chunks_dir=CHUNKS_DIR,
    parse_fn=parse_chunks,
    cells_per_chunk=CELLS_PER_CHUNK,
    total_tbi_cells=TOTAL_START_CELLS,  # field is a misnomer; holds total start-state cells
    goal_state=GOAL_STATE,
)
print(f"Scanning {len(pickles):,} pickles for per-cell, per-gene shifts ...")
per_gene_df = aggregate_per_gene(
    iter_per_cell_shifts(spec),
    donor_per_cell,
    TOTAL_START_CELLS,
    token_to_ensembl,
    token_to_symbol,
)
# aggregate_per_gene hard-codes "percent_TBI_cells_perturbed"; rename for clarity here.
per_gene_df = per_gene_df.rename(
    columns={"percent_TBI_cells_perturbed": f"percent_{START_STATE}_cells_perturbed"}
)
print(f"  Genes with ≥1 perturbation: {len(per_gene_df):,}")

# ── 3. merge with ISP stats (FDR/pval/Sig) ─────────────────────────────────
combined = merge_with_isp_stats(per_gene_df, raw_csv)

# ── 4. emit ranked tables ──────────────────────────────────────────────────
# apply_high_confidence_filter reads "percent_TBI_cells_perturbed" — undo our rename
# just for the filter call, then re-rename to the start-state name.
combined_for_filter = combined.rename(
    columns={f"percent_{START_STATE}_cells_perturbed": "percent_TBI_cells_perturbed"}
)
filtered_raw = apply_high_confidence_filter(
    combined_for_filter, min_cells=MIN_CELLS_PERTURBED, min_pct=MIN_PCT_CELLS
).rename(columns={"percent_TBI_cells_perturbed": f"percent_{START_STATE}_cells_perturbed"})

unfiltered = rank_and_order(combined)
filtered   = rank_and_order(filtered_raw)

unf_path  = STATS_DIR / f"{OUTPUT_PREFIX}_unfiltered_all_gene_ranked.tsv"
filt_path = STATS_DIR / f"{OUTPUT_PREFIX}_filtered_high_confidence_gene_ranked.tsv"
unfiltered.to_csv(unf_path, sep="\t", index=False)
filtered.to_csv(filt_path, sep="\t", index=False)

# ── 5. report ──────────────────────────────────────────────────────────────
print(f"\nOutput files:")
print(f"  {unf_path}  ({len(unfiltered):,} genes)")
print(f"  {filt_path}  ({len(filtered):,} genes, n_cells>={MIN_CELLS_PERTURBED} OR pct>={MIN_PCT_CELLS}%)")

display_cols = [
    "rank", "gene_symbol", "ensembl_id",
    "n_cells_perturbed", f"percent_{START_STATE}_cells_perturbed",
    "mean_shift", "shift_sem",
    "n_donors_contributing", "donor_consistency_score",
    "goal_end_fdr", "sig",
]
display_cols = [c for c in display_cols if c in filtered.columns]

print(f"\nTop 20 candidates by mean_shift toward {GOAL_STATE} (mode={args.mode}, direction={START_STATE}→{GOAL_STATE}):")
if len(filtered) == 0:
    print("  (no genes passed the filter)")
else:
    print(filtered[display_cols].head(20).to_string(index=False))
