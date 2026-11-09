"""
Helpers for computing per-cell, per-gene perturbation statistics and
coverage-based safeguards from Geneformer InSilicoPerturber raw pickles.

Each *_raw.pickle file holds, for one TBI cell, a dict:
  {'TBI':     {(gene_token, 'cell_emb'): [shift_to_TBI_state]},
   'control': {(gene_token, 'cell_emb'): [shift_to_control_state]}}
The shift value is the cosine similarity of (perturbed - original) cell
embedding relative to the corresponding state embedding (verified to match
InSilicoPerturberStats.Shift_to_goal_end when averaged across cells).

This module is shared by the full-run and pilot aggregation scripts.
"""
import pickle
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np
import pandas as pd

GENEFORMER_PKG = Path("/work/pi_cfallini_uri_edu/Geneformer_TBI/models/Geneformer/geneformer")
TOKEN_DICT_PATH = GENEFORMER_PKG / "token_dictionary_gc104M.pkl"
GENE_NAME_DICT_PATH = GENEFORMER_PKG / "gene_name_id_dict_gc104M.pkl"

# Full-run pickle: in_silico_delete_<CELLTYPE>_chunk###_cells####-####_dict_cell_embs_<cell>batch<B>_raw.pickle
FULL_RUN_RE = re.compile(
    r"^in_silico_delete_excitatory_chunk(\d+)_cells\d+-\d+"
    r"_dict_cell_embs_(\d+)batch(-?\d+)_raw\.pickle$"
)
# Pilot pickle: in_silico_delete_pilot_excitatory_nGene1000_dict_cell_embs_<cell>batch<B>_raw.pickle
PILOT_RE = re.compile(
    r"^in_silico_delete_pilot_excitatory_nGene1000"
    r"_dict_cell_embs_(\d+)batch(-?\d+)_raw\.pickle$"
)


def parse_full_run(name: str):
    m = FULL_RUN_RE.match(name)
    if m is None:
        return None
    return int(m.group(1)), int(m.group(2)), int(m.group(3))  # chunk, cell_in_chunk, batch


def parse_pilot(name: str):
    m = PILOT_RE.match(name)
    if m is None:
        return None
    return 0, int(m.group(1)), int(m.group(2))  # chunk=0, cell_idx, batch


def make_parse_full_run(celltype_prefix: str, isp_mode: str = "delete"):
    """Factory: build a parser for chunk-pickle filenames of the form
        in_silico_<isp_mode>_<celltype_prefix>_chunk###_cells####-####
        _dict_cell_embs_<cell>batch<B>_raw.pickle
    isp_mode is "delete" (default, matches the zero-shot pipeline) or "overexpress"
    (used by the fine-tuned control→TBI screen).
    Returns a callable parser(name) -> (chunk_id, cell_in_chunk, batch) or None.
    """
    pat = re.compile(
        rf"^in_silico_{re.escape(isp_mode)}_{re.escape(celltype_prefix)}_chunk(\d+)_cells\d+-\d+"
        rf"_dict_cell_embs_(\d+)batch(-?\d+)_raw\.pickle$"
    )

    def parse(name: str):
        m = pat.match(name)
        if m is None:
            return None
        return int(m.group(1)), int(m.group(2)), int(m.group(3))

    return parse


def load_token_maps() -> tuple[dict[int, str], dict[int, str]]:
    """Returns (token_to_ensembl, token_to_symbol)."""
    with open(TOKEN_DICT_PATH, "rb") as fh:
        token_dict = pickle.load(fh)  # ENSG -> token int
    with open(GENE_NAME_DICT_PATH, "rb") as fh:
        symbol_to_ensembl = pickle.load(fh)  # symbol -> ENSG
    ensembl_to_symbol: dict[str, str] = {}
    for sym, ens in symbol_to_ensembl.items():
        ensembl_to_symbol.setdefault(ens, sym)
    token_to_ensembl = {
        tok: ens for ens, tok in token_dict.items() if isinstance(ens, str) and ens.startswith("ENSG")
    }
    token_to_symbol = {tok: ensembl_to_symbol.get(ens, "") for tok, ens in token_to_ensembl.items()}
    return token_to_ensembl, token_to_symbol


def load_donor_map(tokenized_path: Path, condition: str = "TBI") -> tuple[list[str], list[str]]:
    """Returns (donor_per_cell, sample_per_cell) for cells matching `condition`,
    in the same order InSilicoPerturber iterates them (dataset order after filter).
    Use condition="TBI" for the zero-shot screen, "control" for the fine-tuned
    control→TBI screen.
    """
    from datasets import load_from_disk

    ds = load_from_disk(str(tokenized_path))
    ds_f = ds.filter(lambda r: r["condition"] == condition)
    return list(ds_f["donor"]), list(ds_f["sample"])


def load_tbi_donor_map(tokenized_path: Path) -> tuple[list[str], list[str]]:
    """Backward-compat wrapper: load_donor_map with condition='TBI'."""
    return load_donor_map(tokenized_path, condition="TBI")


@dataclass
class RunSpec:
    chunks_dir: Path
    parse_fn: callable
    cells_per_chunk: int
    total_tbi_cells: int
    goal_state: str = "control"


def iter_per_cell_shifts(spec: RunSpec) -> Iterator[tuple[int, int, float]]:
    """Yield (global_cell_idx, gene_token, shift_to_goal) for every (cell, gene) entry."""
    for pkl in sorted(spec.chunks_dir.glob("*_raw.pickle")):
        parsed = spec.parse_fn(pkl.name)
        if parsed is None:
            continue
        chunk_id, cell_in_chunk, _batch = parsed
        global_idx = chunk_id * spec.cells_per_chunk + cell_in_chunk
        if global_idx >= spec.total_tbi_cells:
            continue
        with open(pkl, "rb") as fh:
            obj = pickle.load(fh)
        inner = obj.get(spec.goal_state, {})
        for (token, kind), vals in inner.items():
            if kind != "cell_emb":
                continue
            for s in vals:
                yield global_idx, int(token), float(s)


def aggregate_per_gene(
    shift_iter: Iterable[tuple[int, int, float]],
    donor_per_cell: list[str],
    total_tbi_cells: int,
    token_to_ensembl: dict[int, str],
    token_to_symbol: dict[int, str],
) -> pd.DataFrame:
    per_gene: dict[int, list[tuple[int, float]]] = defaultdict(list)
    for cell_idx, token, shift in shift_iter:
        per_gene[token].append((cell_idx, shift))

    rows = []
    for token, pairs in per_gene.items():
        cells = np.array([c for c, _ in pairs], dtype=np.int64)
        shifts = np.array([s for _, s in pairs], dtype=np.float64)
        n = shifts.size

        donors = [donor_per_cell[c] for c in cells]
        unique_donors = sorted(set(donors))
        n_donors = len(unique_donors)

        mean_shift = float(shifts.mean())
        if n_donors >= 2 and mean_shift != 0:
            donor_means = {
                d: shifts[np.array([dd == d for dd in donors])].mean() for d in unique_donors
            }
            overall_sign = np.sign(mean_shift)
            sign_match = sum(1 for dm in donor_means.values() if np.sign(dm) == overall_sign)
            donor_consistency = sign_match / n_donors
        else:
            donor_consistency = np.nan

        rows.append(
            {
                "gene_token": int(token),
                "gene_symbol": token_to_symbol.get(token, ""),
                "ensembl_id": token_to_ensembl.get(token, ""),
                "n_cells_perturbed": n,
                "percent_TBI_cells_perturbed": 100.0 * n / total_tbi_cells,
                "mean_shift": mean_shift,
                "median_shift": float(np.median(shifts)),
                "shift_sd": float(shifts.std(ddof=1)) if n > 1 else np.nan,
                "shift_sem": float(shifts.std(ddof=1) / np.sqrt(n)) if n > 1 else np.nan,
                "n_donors_contributing": n_donors,
                "donor_consistency_score": donor_consistency,
            }
        )
    return pd.DataFrame(rows)


def merge_with_isp_stats(per_gene_df: pd.DataFrame, isp_csv: Path) -> pd.DataFrame:
    """Add Goal_end_vs_random_pval / Goal_end_FDR / Sig from InSilicoPerturberStats output."""
    isp = pd.read_csv(isp_csv)
    keep = [c for c in ["Ensembl_ID", "Goal_end_vs_random_pval", "Goal_end_FDR", "Sig"] if c in isp.columns]
    isp = isp[keep].rename(
        columns={
            "Goal_end_vs_random_pval": "goal_end_vs_random_pval",
            "Goal_end_FDR": "goal_end_fdr",
            "Sig": "sig",
        }
    )
    return per_gene_df.merge(isp, left_on="ensembl_id", right_on="Ensembl_ID", how="left").drop(
        columns=["Ensembl_ID"]
    )


OUTPUT_COLS_ORDERED = [
    "rank",
    "gene_symbol",
    "ensembl_id",
    "n_cells_perturbed",
    "percent_TBI_cells_perturbed",
    "mean_shift",
    "median_shift",
    "shift_sd",
    "shift_sem",
    "n_donors_contributing",
    "donor_consistency_score",
    "goal_end_vs_random_pval",
    "goal_end_fdr",
    "sig",
    "gene_token",
]


def rank_and_order(df: pd.DataFrame) -> pd.DataFrame:
    df = df.sort_values("mean_shift", ascending=False).reset_index(drop=True)
    df.insert(0, "rank", range(1, len(df) + 1))
    cols = [c for c in OUTPUT_COLS_ORDERED if c in df.columns]
    extras = [c for c in df.columns if c not in cols]
    return df[cols + extras]


def apply_high_confidence_filter(
    df: pd.DataFrame,
    min_cells: int = 25,
    min_pct: float = 2.0,
    drop_universal_artifacts: bool = True,
) -> pd.DataFrame:
    """High-confidence filter: coverage gate (n_cells_perturbed >= min_cells OR
    percent_TBI_cells_perturbed >= min_pct), optionally minus universal abundance
    artifacts (MT-* mitochondrial-encoded, RPL*/RPS* cytoplasmic ribosomal,
    MALAT1, hemoglobins).

    Cell-type-specific contamination markers (neuronal/astrocyte/oligo/microglia
    cross-contamination) are NOT removed here — those belong to the biological-
    candidate filter, where the choice of contamination set depends on the
    target cell type.
    """
    mask = (df["n_cells_perturbed"] >= min_cells) | (df["percent_TBI_cells_perturbed"] >= min_pct)
    result = df[mask].copy()
    if drop_universal_artifacts:
        from biological_gene_classes import ARTIFACT_LITERALS, HEMOGLOBINS

        sym = result["gene_symbol"].fillna("").astype(str)
        keep = ~(
            sym.str.startswith("MT-")
            | sym.str.startswith("RPL")
            | sym.str.startswith("RPS")
            | sym.isin(ARTIFACT_LITERALS)
            | sym.isin(HEMOGLOBINS)
        )
        result = result[keep].copy()
    return result


def compute_detection_fractions(
    tokenized_path,
    state_key: str,
    states: list[str],
    cell_type_key: str | None = "cell_type",
):
    """Per-gene detection fraction = fraction of cells whose rank-value encoding
    (input_ids) contains the gene token, computed within each condition `state`
    and (if cell_type_key given) within each (cell_type, state).

    This is the metric behind the Geneformer preprint eligibility rule
    (2024.08.16.608180v1, Methods): "tested in silico perturbation of each gene
    that was detected in at least 20% of the start cells or goal end cells."
    Only the requested `states` are scanned (pass [start_state, goal_state]).

    Returns (det_pooled, det_byct):
      det_pooled[state][token]            -> fraction in cells where state_key==state
      det_byct[(cell_type, state)][token] -> fraction within that cell type
    det_byct is {} when cell_type_key is None.
    """
    from collections import Counter, defaultdict
    from datasets import load_from_disk

    ds = load_from_disk(str(tokenized_path))
    state_vals = list(ds[state_key])
    ct_vals = list(ds[cell_type_key]) if cell_type_key else None
    input_ids = ds["input_ids"]
    want = set(states)

    pooled_cnt = {s: Counter() for s in states}
    pooled_n = {s: 0 for s in states}
    byct_cnt: dict = defaultdict(Counter)
    byct_n: dict = defaultdict(int)

    for i, st in enumerate(state_vals):
        if st not in want:
            continue
        toks = set(input_ids[i])
        pooled_cnt[st].update(toks)
        pooled_n[st] += 1
        if ct_vals is not None:
            key = (ct_vals[i], st)
            byct_cnt[key].update(toks)
            byct_n[key] += 1

    det_pooled = {
        s: ({t: c / pooled_n[s] for t, c in pooled_cnt[s].items()} if pooled_n[s] else {})
        for s in states
    }
    det_byct = {
        key: {t: c / byct_n[key] for t, c in cnt.items()}
        for key, cnt in byct_cnt.items()
    }
    return det_pooled, det_byct


def apply_eligibility_filter(
    df: pd.DataFrame,
    det_start: dict,
    det_goal: dict,
    thresh: float = 0.20,
    drop_universal_artifacts: bool = True,
):
    """Geneformer-preprint gene eligibility filter: keep genes detected in
    >= `thresh` of start OR goal cells, optionally minus universal abundance
    artifacts (MT-*, RPL*/RPS*, MALAT1, hemoglobins). Annotates det_start /
    det_goal / eligible_pct onto a COPY of every input row.

    Returns (annotated_all, filtered):
      annotated_all carries the det columns for ALL genes (use for the
        unfiltered table);
      filtered is the eligible, non-artifact subset (the high-confidence table).
    """
    annotated = df.copy()
    annotated["det_start"] = annotated["gene_token"].map(lambda t: det_start.get(int(t), 0.0))
    annotated["det_goal"] = annotated["gene_token"].map(lambda t: det_goal.get(int(t), 0.0))
    annotated["eligible_pct"] = (
        (annotated["det_start"] >= thresh) | (annotated["det_goal"] >= thresh)
    )

    filtered = annotated[annotated["eligible_pct"]].copy()
    if drop_universal_artifacts:
        from biological_gene_classes import ARTIFACT_LITERALS, HEMOGLOBINS

        sym = filtered["gene_symbol"].fillna("").astype(str)
        keep = ~(
            sym.str.startswith("MT-")
            | sym.str.startswith("RPL")
            | sym.str.startswith("RPS")
            | sym.isin(ARTIFACT_LITERALS)
            | sym.isin(HEMOGLOBINS)
        )
        filtered = filtered[keep].copy()
    return annotated, filtered
