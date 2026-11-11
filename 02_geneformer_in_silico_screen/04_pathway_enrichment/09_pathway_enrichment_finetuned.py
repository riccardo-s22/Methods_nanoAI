#!/usr/bin/env python
"""Pathway enrichment (ORA) for the fine-tuned control→TBI screen.

Runs `gseapy.enrichr` against 4 libraries with the screen's tested-gene
universe as background, for two candidate sets per mode:
  • TIER_HIGH    — high-priority candidates (gene_class != "other",
                    FDR<1e-10, n_cells>=50, consist>=0.85)
  • TIER_ALL_BIO — all biological-filter-passing candidates (broader view
                    including "other" class genes)

Outputs:
  <run_root>/enrichment_<mode>/
      high_tier/<library>.tsv
      high_tier/input_gene_list.txt
      all_bio/<library>.tsv
      all_bio/input_gene_list.txt
      summary_top5.tsv

Per-celltype config (run_root, celltype_prefix) lives in
`finetuned_celltype_config.py`. Default celltype is `excit`.

Usage:
  python scripts/09_pathway_enrichment_finetuned.py --mode delete      [--celltype excit|oligo]
  python scripts/09_pathway_enrichment_finetuned.py --mode overexpress [--celltype excit|oligo]

In this repository the script is 04_pathway_enrichment/09_pathway_enrichment_finetuned.py
and <run_root> is $GF_WORK/results/perturb_neurons_finetuned_<direction> (see
../gf_paths.py). It needs internet access (Enrichr).
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import gseapy as gp
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
# Section path configuration: see ../gf_paths.py. finetuned_celltype_config.py
# sits in 02_model_finetuning/ft_excitatory in this repository (the scripts were
# written for one flat scripts/ folder); add_ft_code_to_syspath() makes it importable.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gf_paths  # noqa: E402

gf_paths.add_ft_code_to_syspath()
from finetuned_celltype_config import CELLTYPE_CONFIG, get_config  # noqa: E402

LIBRARIES = [
    "GO_Biological_Process_2023",
    "Reactome_2022",
    "KEGG_2021_Human",
    "MSigDB_Hallmark_2020",
]


def overlap_count(row: pd.Series) -> int:
    if "Overlap" in row:
        return row["Overlap"]
    genes = row.get("Genes", "")
    return len([g for g in str(genes).split(";") if g])


def run_enrichr(gene_list: list[str], bg_genes: list[str], tag_dir: Path,
                summary_rows: list[dict], tier_label: str) -> None:
    tag_dir.mkdir(parents=True, exist_ok=True)
    (tag_dir / "input_gene_list.txt").write_text("\n".join(gene_list) + "\n")

    for lib in LIBRARIES:
        print(f"  → {lib}", flush=True)
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=lib,
                organism="human",
                background=bg_genes,
                outdir=None,
                cutoff=1.0,
            )
        except Exception as exc:
            print(f"    !! {lib} failed: {exc}")
            continue

        df = enr.results.copy().sort_values("Adjusted P-value")
        out_path = tag_dir / f"{lib}.tsv"
        df.to_csv(out_path, sep="\t", index=False)
        n_sig = int((df["Adjusted P-value"] < 0.05).sum())
        print(f"    wrote {out_path.name}: {len(df):,} terms, {n_sig} with FDR < 0.05")

        for _, row in df[df["Adjusted P-value"] < 0.05].head(5).iterrows():
            summary_rows.append({
                "tier":         tier_label,
                "library":      lib,
                "term":         row["Term"],
                "overlap":      overlap_count(row),
                "adjusted_p":   row["Adjusted P-value"],
                "odds_ratio":   row.get("Odds Ratio", float("nan")),
                "combined":     row.get("Combined Score", float("nan")),
                "genes":        row.get("Genes", ""),
            })


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["delete", "overexpress"], required=True)
    ap.add_argument("--celltype", choices=sorted(CELLTYPE_CONFIG.keys()),
                    default="excit",
                    help="Cell type config from finetuned_celltype_config.py (default: excit)")
    args = ap.parse_args()

    cfg = get_config(args.celltype)
    RUN_ROOT = cfg["run_root"]
    CELLTYPE_PREFIX = cfg["celltype_prefix"]

    STATS_DIR = RUN_ROOT / f"stats_{args.mode}"
    CAND_TSV  = STATS_DIR / f"{CELLTYPE_PREFIX}_{args.mode}_filtered_biological_candidate_ranked.tsv"
    UNIV_TSV  = STATS_DIR / f"{CELLTYPE_PREFIX}_{args.mode}_unfiltered_all_gene_ranked.tsv"
    OUT_DIR   = RUN_ROOT / f"enrichment_{args.mode}"

    if not CAND_TSV.exists():
        raise SystemExit(f"ERROR: candidates TSV not found: {CAND_TSV}\n"
                         f"Run 07_biological_candidate_filter_finetuned.py --mode {args.mode} --celltype {args.celltype} first.")
    if not UNIV_TSV.exists():
        raise SystemExit(f"ERROR: universe TSV not found: {UNIV_TSV}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"=== Pathway enrichment | celltype={args.celltype} mode={args.mode} ===")
    print(f"Candidates: {CAND_TSV.name}")
    print(f"Universe:   {UNIV_TSV.name}")
    print(f"Outdir:     {OUT_DIR}")

    cands = pd.read_csv(CAND_TSV, sep="\t")
    universe = pd.read_csv(UNIV_TSV, sep="\t")
    bg_genes = sorted({g for g in universe["gene_symbol"].dropna() if g})

    cands["gene_symbol"] = cands["gene_symbol"].fillna("").astype(str)
    high = cands[cands["biological_priority"] == "high"].copy()
    high_genes = sorted({g for g in high["gene_symbol"].dropna() if g})
    all_genes  = sorted({g for g in cands["gene_symbol"].dropna() if g})
    print(f"\nUniverse (tested genes): {len(bg_genes):,}")
    print(f"High-tier candidates:    {len(high_genes):,}")
    print(f"All biological cand.:    {len(all_genes):,}")

    summary_rows: list[dict] = []

    print(f"\n--- High-tier ({len(high_genes)} genes) ---")
    run_enrichr(high_genes, bg_genes, OUT_DIR / "high_tier", summary_rows, "high")

    print(f"\n--- All biological ({len(all_genes)} genes) ---")
    run_enrichr(all_genes, bg_genes, OUT_DIR / "all_bio", summary_rows, "all_bio")

    summary = pd.DataFrame(summary_rows)
    summary_path = OUT_DIR / "summary_top5.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"\nWrote summary: {summary_path}")
    print(f"  high-tier sig terms (top 5/lib): "
          f"{int((summary['tier']=='high').sum())}")
    print(f"  all-bio sig terms (top 5/lib):  "
          f"{int((summary['tier']=='all_bio').sum())}")


if __name__ == "__main__":
    main()
