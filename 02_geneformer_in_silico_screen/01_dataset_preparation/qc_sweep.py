#!/usr/bin/env python
"""Sweep several MAD/QC configurations against the paper's 15,110-nuclei target,
to identify the configuration that best reproduces the paper's filtering."""
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation

BASE = Path("/mnt/d/TBI_SC/Geneformer_TBI")
SRC = BASE / "h5ad/Garza_GSE209552_all_raw_counts_geneformer_ready.h5ad"
ANN = BASE / "h5ad/Garza_GSE209552_all_raw_counts_paperQC_annotated.h5ad"
PAPER_TOTAL = 15110; PAPER_TBI = 6806; PAPER_CONTROL = 8304

adata = sc.read_h5ad(SRC)
ann = sc.read_h5ad(ANN)
adata.obs["cell_type_prev"] = pd.Series(ann.obs["cell_type"].astype(str).values,
                                        index=ann.obs_names).reindex(adata.obs_names).fillna("removed_qc").values
print(f"Input cells: {adata.n_obs}\n")

def mad_filter_per_sample(adata, nmads_mt, nmads_g, log_g):
    keep = np.zeros(adata.n_obs, dtype=bool)
    for s, idx in adata.obs.groupby("sample", observed=False).groups.items():
        pct_mt = adata.obs.loc[idx, "pct_counts_mt"].astype(float).values
        n_g = adata.obs.loc[idx, "n_genes_by_counts"].astype(float).values
        g_metric = np.log10(np.clip(n_g, 1, None)) if log_g else n_g
        med_mt = np.median(pct_mt); mad_mt = median_abs_deviation(pct_mt) or 1e-6
        med_g  = np.median(g_metric); mad_g  = median_abs_deviation(g_metric) or 1e-6
        sel = (pct_mt <= med_mt + nmads_mt * mad_mt) & (g_metric >= med_g - nmads_g * mad_g)
        keep[adata.obs_names.get_indexer(idx)] = sel
    return keep

def mad_filter_global(adata, nmads_mt, nmads_g, log_g):
    pct_mt = adata.obs["pct_counts_mt"].astype(float).values
    n_g = adata.obs["n_genes_by_counts"].astype(float).values
    g_metric = np.log10(np.clip(n_g, 1, None)) if log_g else n_g
    med_mt = np.median(pct_mt); mad_mt = median_abs_deviation(pct_mt) or 1e-6
    med_g  = np.median(g_metric); mad_g  = median_abs_deviation(g_metric) or 1e-6
    return (pct_mt <= med_mt + nmads_mt * mad_mt) & (g_metric >= med_g - nmads_g * mad_g)

n_g_arr = adata.obs["n_genes_by_counts"].astype(float).values
cond = adata.obs["condition"].astype(str).values

configs = [
    ("nGene>=1000 only (no MAD)",            np.ones(adata.n_obs, dtype=bool)),
    ("per-sample MAD-3 lin + nGene>=1000",   mad_filter_per_sample(adata, 3, 3, log_g=False)),
    ("per-sample MAD-3 log + nGene>=1000",   mad_filter_per_sample(adata, 3, 3, log_g=True)),
    ("per-sample MAD-5 log + nGene>=1000",   mad_filter_per_sample(adata, 5, 5, log_g=True)),
    ("global MAD-3 log + nGene>=1000",       mad_filter_global(adata, 3, 3, log_g=True)),
    ("global MAD-5 log + nGene>=1000",       mad_filter_global(adata, 5, 5, log_g=True)),
]

print(f"{'config':50s} {'total':>6s}  {'TBI':>6s}  {'ctrl':>6s}  vs paper")
print(f"{'paper':50s} {PAPER_TOTAL:>6d}  {PAPER_TBI:>6d}  {PAPER_CONTROL:>6d}")
print("-" * 90)
rows = []
for label, mad_mask in configs:
    final = mad_mask & (n_g_arr >= 1000)
    n = int(final.sum())
    n_tbi = int(((cond == "TBI") & final).sum())
    n_ctrl = int(((cond == "control") & final).sum())
    delta = n - PAPER_TOTAL
    print(f"{label:50s} {n:>6d}  {n_tbi:>6d}  {n_ctrl:>6d}  Δ={delta:+d}")
    rows.append({"config": label, "total": n, "TBI": n_tbi, "control": n_ctrl, "delta_vs_paper": delta})

out = BASE / "reports/qc_sweep.tsv"
pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
print(f"\nWrote {out}")
