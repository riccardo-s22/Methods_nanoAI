#!/usr/bin/env python
"""
Biological-candidate filter for the fine-tuned control→TBI Geneformer
perturbation results (both delete and overexpress modes).

Mirrors `07_biological_candidate_filter.py` (zero-shot excitatory primary
screen) but reads the fine-tuned stats dirs and supports both perturb
modes via --mode.

Source (per --mode + --celltype):
  <run_root>/stats_<mode>/<celltype_prefix>_<mode>_unfiltered_all_gene_ranked.tsv

Outputs (same stats_<mode> dir):
  <celltype_prefix>_<mode>_filtered_biological_candidate_ranked.tsv
  <celltype_prefix>_<mode>_top100_biological_candidate_summary.tsv

Per-celltype config (run_root, celltype_prefix) lives in
`finetuned_celltype_config.py`. Default celltype is `excit`.

Filters (identical to the primary screen — biology-driven, mode-agnostic):
  1. mean_shift > 0   (perturbation moves control toward TBI)
  2. n_cells_perturbed >= 25
  3. n_donors_contributing >= 3
  4. donor_consistency_score >= 0.7
  5. drop MT-       (mitochondrial-encoded)
  6. drop RPL/RPS   (cytoplasmic ribosomal proteins)
  7. drop literals  (MALAT1)
  8. drop hemoglobins
  9. drop non-neuronal contamination markers
 10. keep IEGs (FOS/JUN/JUNB/ATF3) if they pass; annotate as stress.

Interpretation notes:
  • delete mode      → loss-of-function TBI drivers (gene loss INDUCES TBI)
  • overexpress mode → gain-of-function TBI drivers (gene OE INDUCES TBI)

Usage:
  python scripts/07_biological_candidate_filter_finetuned.py --mode delete      [--celltype excit|oligo]
  python scripts/07_biological_candidate_filter_finetuned.py --mode overexpress [--celltype excit|oligo]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from biological_gene_classes import (
    ARTIFACT_LITERALS,
    HEMOGLOBINS,
    NONNEURONAL_MARKERS,
    classify_gene,
)
from finetuned_celltype_config import CELLTYPE_CONFIG, get_config

MIN_N_CELLS       = 25
MIN_N_DONORS      = 3
MIN_DONOR_CONSIST = 0.7

IEG_KEEPLIST = {"FOS", "JUN", "JUNB", "ATF3"}

INTERPRETATION = {
    "delete":      "Top hits: gene loss INDUCES TBI state (loss-of-function TBI drivers)",
    "overexpress": "Top hits: gene over-expression INDUCES TBI state (gain-of-function TBI drivers)",
}


def filter_step(df: pd.DataFrame, mask: pd.Series, label: str) -> pd.DataFrame:
    kept = int(mask.sum())
    dropped = len(df) - kept
    print(f"  {label}: keep {kept:,} / drop {dropped:,}  (remaining: {kept:,})")
    return df[mask].copy()


def assign_artifact_flag(symbol: str) -> str:
    """Only IEG-kept genes are annotated. Other artifacts are removed
    before reaching this function."""
    if symbol in IEG_KEEPLIST:
        return "IEG_stress_response"
    return ""


def assign_biological_priority(row: pd.Series) -> str:
    """Tier logic recalibrated for the fine-tuned control→TBI direction.

    The primary-screen tiers (`07_biological_candidate_filter.py`) required
    n_donors >= 4/6, sized for 12 GARZA TBI start-state donors. This screen
    starts from CONTROL cells with only 3 donors total, so the donor count
    is capped at 3 and that requirement is unreachable. We drop the
    n_donors gate (the prefilters already require n_donors >= 3 AND
    donor_consistency >= 0.7, which for 3 donors means all-3-agree) and
    discriminate tiers by signal strength (FDR + n_cells) and biological
    actionability (gene_class).
    """
    fdr = row.get("goal_end_fdr", np.nan)
    consist = row.get("donor_consistency_score", np.nan)
    n_cells = row.get("n_cells_perturbed", 0)
    gene_class = row.get("gene_class", "other")

    fdr = float("inf") if pd.isna(fdr) else float(fdr)
    consist = 0.0 if pd.isna(consist) else float(consist)

    if (fdr < 1e-10 and consist >= 0.85 and n_cells >= 50
            and gene_class != "other"):
        return "high"
    if (fdr < 1e-5 and consist >= 0.75
            and (n_cells >= 100 or gene_class != "other")):
        return "medium"
    return "low"


def assign_reason(row: pd.Series, pct_col: str) -> str:
    parts = []
    if row["artifact_flag"] == "IEG_stress_response":
        parts.append("IEG/stress-response (kept per spec)")
    cls = row["gene_class"]
    if cls != "other":
        parts.append(f"class={cls}")
    parts.append(
        f"n_cells={int(row['n_cells_perturbed'])}, "
        f"donors={int(row['n_donors_contributing'])}, "
        f"consistency={row['donor_consistency_score']:.2f}"
    )
    fdr = row.get("goal_end_fdr", np.nan)
    if pd.notna(fdr):
        parts.append(f"FDR={fdr:.1e}")
    return "; ".join(parts)


# ── CLI ───────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--mode", choices=["delete", "overexpress"], required=True)
parser.add_argument("--celltype", choices=sorted(CELLTYPE_CONFIG.keys()),
                    default="excit",
                    help="Cell type config from finetuned_celltype_config.py (default: excit)")
args = parser.parse_args()

cfg = get_config(args.celltype)
RUN_ROOT = cfg["run_root"]
CELLTYPE_PREFIX = cfg["celltype_prefix"]

STATS_DIR = RUN_ROOT / f"stats_{args.mode}"
OUTPUT_PREFIX = f"{CELLTYPE_PREFIX}_{args.mode}"
SOURCE   = STATS_DIR / f"{OUTPUT_PREFIX}_unfiltered_all_gene_ranked.tsv"
OUT_RANK = STATS_DIR / f"{OUTPUT_PREFIX}_filtered_biological_candidate_ranked.tsv"
OUT_TOP  = STATS_DIR / f"{OUTPUT_PREFIX}_top100_biological_candidate_summary.tsv"

print(f"=== Fine-tuned control→TBI biological filter | celltype={args.celltype} mode={args.mode} ===")
print(f"{INTERPRETATION[args.mode]}")

# ── 1. load source ────────────────────────────────────────────────────────
if not SOURCE.exists():
    sys.exit(
        f"ERROR: source TSV not found: {SOURCE}\n"
        f"Run scripts/05_get_perturb_stats_finetuned.py --mode {args.mode} first."
    )

df = pd.read_csv(SOURCE, sep="\t")
print(f"\nLoaded {len(df):,} genes from {SOURCE.name}")

# The fine-tuned aggregator renames percent_TBI_cells_perturbed →
# percent_control_cells_perturbed (since start=control). Accept either.
if "percent_control_cells_perturbed" in df.columns:
    pct_col = "percent_control_cells_perturbed"
elif "percent_TBI_cells_perturbed" in df.columns:
    pct_col = "percent_TBI_cells_perturbed"
else:
    sys.exit("ERROR: source TSV missing percent_*_cells_perturbed column")

required = {"gene_symbol", "ensembl_id", "n_cells_perturbed",
            pct_col, "mean_shift", "median_shift",
            "shift_sd", "shift_sem", "n_donors_contributing",
            "donor_consistency_score", "goal_end_fdr", "sig"}
missing = required - set(df.columns)
if missing:
    sys.exit(f"ERROR: source TSV missing required columns: {sorted(missing)}")

df["gene_symbol"] = df["gene_symbol"].fillna("").astype(str)
df["donor_consistency_score"] = pd.to_numeric(df["donor_consistency_score"], errors="coerce")

# ── 2. apply filters ──────────────────────────────────────────────────────
print("\nApplying biological filters:")

df = filter_step(df, df["mean_shift"] > 0, "mean_shift > 0  (control→TBI direction)")
df = filter_step(df, df["n_cells_perturbed"] >= MIN_N_CELLS,
                 f"n_cells_perturbed >= {MIN_N_CELLS}")
df = filter_step(df, df["n_donors_contributing"] >= MIN_N_DONORS,
                 f"n_donors_contributing >= {MIN_N_DONORS}")
df = filter_step(df, df["donor_consistency_score"] >= MIN_DONOR_CONSIST,
                 f"donor_consistency_score >= {MIN_DONOR_CONSIST}")

mt_mask = ~df["gene_symbol"].str.startswith("MT-", na=False)
df = filter_step(df, mt_mask, "drop MT- mitochondrial-encoded")

rib_mask = ~(df["gene_symbol"].str.startswith("RPL", na=False)
             | df["gene_symbol"].str.startswith("RPS", na=False))
df = filter_step(df, rib_mask, "drop RPL/RPS ribosomal")

lit_mask = ~df["gene_symbol"].isin(ARTIFACT_LITERALS)
df = filter_step(df, lit_mask, f"drop {sorted(ARTIFACT_LITERALS)}")

hb_mask = ~df["gene_symbol"].isin(HEMOGLOBINS)
df = filter_step(df, hb_mask, f"drop hemoglobins ({len(HEMOGLOBINS)} symbols)")

nn_mask = ~df["gene_symbol"].isin(NONNEURONAL_MARKERS)
df = filter_step(df, nn_mask, f"drop non-neuronal markers ({len(NONNEURONAL_MARKERS)} symbols)")

# ── 3. annotate ───────────────────────────────────────────────────────────
df["artifact_flag"] = df["gene_symbol"].map(assign_artifact_flag)
df["gene_class"] = df["gene_symbol"].map(classify_gene)
df["biological_priority"] = df.apply(assign_biological_priority, axis=1)
df["reason"] = df.apply(lambda r: assign_reason(r, pct_col), axis=1)

for ieg in sorted(IEG_KEEPLIST):
    if ieg in set(df["gene_symbol"]):
        print(f"  IEG kept: {ieg} -> flagged 'IEG_stress_response'")

# ── 4. re-rank and write ──────────────────────────────────────────────────
df = df.sort_values("mean_shift", ascending=False).reset_index(drop=True)
df["rank"] = range(1, len(df) + 1)

ordered = [
    "rank", "gene_symbol", "ensembl_id",
    "artifact_flag", "gene_class", "biological_priority", "reason",
    "n_cells_perturbed", pct_col,
    "mean_shift", "median_shift", "shift_sd", "shift_sem",
    "n_donors_contributing", "donor_consistency_score",
    "goal_end_vs_random_pval", "goal_end_fdr", "sig", "gene_token",
]
ordered = [c for c in ordered if c in df.columns]
extras = [c for c in df.columns if c not in ordered]
df = df[ordered + extras]

df.to_csv(OUT_RANK, sep="\t", index=False)
df.head(100).to_csv(OUT_TOP, sep="\t", index=False)

print(f"\nOutput files:")
print(f"  {OUT_RANK}  ({len(df):,} genes)")
print(f"  {OUT_TOP}  ({min(100, len(df)):,} genes)")

# ── 5. distribution summaries ─────────────────────────────────────────────
class_counts = df["gene_class"].value_counts()
print(f"\nGene-class distribution ({args.mode} biological candidates):")
for cls, n in class_counts.items():
    print(f"  {cls:32s} {n:5d}")

prio_counts = df["biological_priority"].value_counts()
print("\nBiological-priority distribution:")
for prio, n in prio_counts.items():
    print(f"  {prio:32s} {n:5d}")

# ── 6. top 30 ─────────────────────────────────────────────────────────────
print(f"\nTop 30 biological candidate genes ({args.mode} mode):")
display_cols = [
    "rank", "gene_symbol", "mean_shift", "n_cells_perturbed",
    "n_donors_contributing", "donor_consistency_score",
    "goal_end_fdr", "gene_class", "reason",
]
display_cols = [c for c in display_cols if c in df.columns]

with pd.option_context("display.max_colwidth", 80, "display.width", 240):
    print(df[display_cols].head(30).to_string(index=False))
