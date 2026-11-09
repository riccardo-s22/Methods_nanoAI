#!/usr/bin/env python
"""Refine the QC pipeline to better match the paper's 15,110 nuclei.

Changes vs the previous (paperQC + scrublet + nGene>=1000) build:
  - DROP scrublet (the paper does not run a doublet caller; their R pipeline
    uses a rat-spike-in contamination filter that does not apply to us).
  - MAD applied to log10(n_genes_by_counts) lower tail (this matches the
    scater::isOutlier default behaviour for skewed count data, which is what
    the paper's pipeline ultimately calls into).
  - MAD applied to pct_counts_mt upper tail in linear space.
  - DROP the n_genes >= 200 floor (paper has FILTER_MIN_NUM_GENES_EXPRESSED=0).
  - Then apply the explicit n_genes >= 1000 hard cutoff the paper documents.
  - Reuse the existing cell_type / leiden annotation produced on the QC-pass
    cells (it is preserved on the obs index of the merged source matrix).

Writes a sweep table for transparency (which step gates how many cells), then
saves the chosen primary object + cell-type subsets.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import median_abs_deviation

BASE = Path("/mnt/d/TBI_SC/Geneformer_TBI")
SRC = BASE / "h5ad/Garza_GSE209552_all_raw_counts_geneformer_ready.h5ad"   # full merged
ANN = BASE / "h5ad/Garza_GSE209552_all_raw_counts_paperQC_annotated.h5ad"  # has cell_type
OUT = BASE / "h5ad/Garza_GSE209552_all_raw_counts_paperQC_nGene1000_annotated.h5ad"
SUBSET_DIR = BASE / "h5ad/celltype_subsets_paperQC_nGene1000"
REPORT = BASE / "reports/paperQC_nGene1000_report.txt"

PAPER_TOTAL = 15110; PAPER_TBI = 6806; PAPER_CONTROL = 8304
NMADS = 3.0

print(f"Reading source merged matrix: {SRC.name}", flush=True)
adata = sc.read_h5ad(SRC)
n0 = adata.n_obs
print(f"  Input: {n0} cells x {adata.n_vars} genes", flush=True)

# Ensure QC metrics.
if "pct_counts_mt" not in adata.obs.columns or "n_genes_by_counts" not in adata.obs.columns:
    adata.var["mt"] = adata.var["gene_symbol"].astype(str).str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)

# ---- Per-sample scater::isOutlier-style MAD (log on n_genes, linear on pct_mt) ----
print("\nPer-sample MAD QC: 3 MAD log10(n_genes) lower, 3 MAD pct_mt upper", flush=True)
keep_mad = np.zeros(adata.n_obs, dtype=bool)
keep_mt  = np.zeros(adata.n_obs, dtype=bool)
keep_g   = np.zeros(adata.n_obs, dtype=bool)
for s, idx in adata.obs.groupby("sample", observed=False).groups.items():
    pct_mt = adata.obs.loc[idx, "pct_counts_mt"].astype(float).values
    n_g = adata.obs.loc[idx, "n_genes_by_counts"].astype(float).values
    log_g = np.log10(np.clip(n_g, 1, None))
    med_mt = np.median(pct_mt); mad_mt = median_abs_deviation(pct_mt) or 1e-6
    med_lg = np.median(log_g);  mad_lg = median_abs_deviation(log_g) or 1e-6
    pass_mt = pct_mt <= med_mt + NMADS * mad_mt
    pass_g  = log_g >= med_lg - NMADS * mad_lg
    sel = pass_mt & pass_g
    positions = adata.obs_names.get_indexer(idx)
    keep_mad[positions] = sel
    keep_mt[positions]  = pass_mt
    keep_g[positions]   = pass_g

adata.obs["pass_mt_mad3"]    = keep_mt.astype(int)
adata.obs["pass_log_g_mad3"] = keep_g.astype(int)
adata.obs["pass_mad_combo"]  = keep_mad.astype(int)

n_mad = int(keep_mad.sum())
print(f"  pass pct_mt MAD only:        {int(keep_mt.sum())}", flush=True)
print(f"  pass log10(n_g) MAD only:    {int(keep_g.sum())}", flush=True)
print(f"  pass combined MAD:           {n_mad}", flush=True)

# Hard floor n_genes >= 1000.
n_g_arr = adata.obs["n_genes_by_counts"].astype(float).values
pass_1000 = n_g_arr >= 1000
adata.obs["pass_nGene_1000"] = pass_1000.astype(int)
final_keep = keep_mad & pass_1000
print(f"  pass n_genes >= 1000:        {int(pass_1000.sum())}", flush=True)
print(f"  pass MAD + n_genes>=1000:    {int(final_keep.sum())}", flush=True)
adata.obs["qc_pass_final_v2"] = final_keep.astype(int)

# ---- Pull in the existing cell_type annotation from the previous build ----
print("\nMerging cell_type from previous paper-QC annotated object ...", flush=True)
ann = sc.read_h5ad(ANN)
ct_map = pd.Series(ann.obs["cell_type"].astype(str).values, index=ann.obs_names)
adata.obs["cell_type"] = ct_map.reindex(adata.obs_names).fillna("removed_qc").values

# Final = MAD pass + nGene>=1000.
final = adata[final_keep].copy()
final.obs["n_counts"] = np.asarray(final.X.sum(axis=1)).ravel().astype(np.float32)
final.obs["filter_pass"] = 1

# Any cells that pass the new QC but were not in the previous annotation set
# (because they were removed by scrublet last round) will have cell_type =
# "removed_qc". We can re-annotate those quickly via nearest-cluster centroid,
# but a simpler and more honest approach is to label them "unannotated_newly_pass"
# and exclude them from cell-type subsets.
unann = (final.obs["cell_type"] == "removed_qc").sum()
print(f"  cells passing new QC but lacking prior annotation: {int(unann)}", flush=True)
final.obs.loc[final.obs["cell_type"] == "removed_qc", "cell_type"] = "unannotated"

assert sparse.issparse(final.X)
x = final.X.data[:200000]
assert np.allclose(x, np.round(x)), "Raw counts must be integer-like."

final.write_h5ad(OUT, compression="gzip")
print(f"\nWrote {OUT}  ({final.n_obs} cells)", flush=True)

# ---------------------------------------------------------------------
# Cell-type subsets.
SUBSET_DIR.mkdir(parents=True, exist_ok=True)
# wipe stale subset files but keep the directory
for p in SUBSET_DIR.glob("*.h5ad"):
    p.unlink()
SUBSET_MAP = {
    "excitatory_neurons_raw_geneformer_ready.h5ad": ["excitatory_neuron"],
    "inhibitory_neurons_raw_geneformer_ready.h5ad": ["inhibitory_neuron"],
    "oligodendrocyte_raw_geneformer_ready.h5ad":    ["oligodendrocyte"],
    "OPC_raw_geneformer_ready.h5ad":                ["OPC"],
}
subset_info = {}
for fname, labels in SUBSET_MAP.items():
    sub = final[final.obs["cell_type"].astype(str).isin(labels)].copy()
    sub.obs["n_counts"] = np.asarray(sub.X.sum(axis=1)).ravel().astype(np.float32)
    sub.obs["filter_pass"] = 1
    sub.write_h5ad(SUBSET_DIR / fname, compression="gzip")
    subset_info[fname] = (sub.n_obs, sub.obs["condition"].value_counts(dropna=False).to_dict())

# ---------------------------------------------------------------------
# Excitatory-neuron gene coverage and Geneformer V2 vocab overlap.
import pickle
V2_PKL = Path("/home/riccardo-s/miniconda3/envs/geneformer/lib/python3.10/site-packages/geneformer/token_dictionary_gc104M.pkl")
exc = final[final.obs["cell_type"] == "excitatory_neuron"]
exc_present = np.asarray((exc.X > 0).sum(axis=0)).ravel()
exc_detected_ens = final.var["ensembl_id"].astype(str).values[exc_present >= 1]
n_exc_detected = int(len(exc_detected_ens))
with open(V2_PKL, "rb") as fh:
    td = pickle.load(fh)
v2_ens = {k for k in td if isinstance(k, str) and k.startswith("ENSG")}
n_in_vocab = len(set(exc_detected_ens) & v2_ens)

# ---------------------------------------------------------------------
# Report.
lines = []
lines.append("=" * 72)
lines.append("Paper-matched primary object (refined): MAD(log10 n_genes, pct_mt) + n>=1000")
lines.append("=" * 72)
lines.append(f"Input cells:                  {n0}")
lines.append(f"After MAD-3 (log10 n_genes,pct_mt): {n_mad}")
lines.append(f"After + n_genes >= 1000:      {int(final_keep.sum())}  "
             f"(paper: {PAPER_TOTAL})")
tbi_n = int((final.obs['condition'] == 'TBI').sum())
ctrl_n = int((final.obs['condition'] == 'control').sum())
lines.append(f"  TBI:                        {tbi_n}  (paper: {PAPER_TBI})")
lines.append(f"  control:                    {ctrl_n}  (paper: {PAPER_CONTROL})")
lines.append("")
lines.append("Cells per sample:")
lines.append(final.obs["sample"].value_counts().sort_index().to_string())
lines.append("")
lines.append("Cells per cell_type:")
lines.append(final.obs["cell_type"].value_counts(dropna=False).to_string())
lines.append("")
lines.append("Cells per (cell_type x condition):")
lines.append(pd.crosstab(final.obs["cell_type"], final.obs["condition"]).to_string())
lines.append("")
lines.append("Core cell-type counts:")
for ct in ["excitatory_neuron", "inhibitory_neuron", "oligodendrocyte", "OPC", "astrocyte", "microglia_macrophage"]:
    n = int((final.obs["cell_type"] == ct).sum())
    lines.append(f"  {ct:24s} {n}")
lines.append("")
lines.append("Excitatory-neuron gene coverage:")
lines.append(f"  Total genes in matrix:               {final.n_vars}")
lines.append(f"  Genes detected (>=1 read >=1 cell):  {n_exc_detected}")
lines.append(f"  Geneformer V2 vocab size (ENSG):     {len(v2_ens)}")
lines.append(f"  Excitatory-detected in V2 vocab:     {n_in_vocab}")
lines.append("")
lines.append("Subset h5ad outputs:")
for fname, (n, cc) in subset_info.items():
    lines.append(f"  {fname}: n={n}, conditions={cc}")
text = "\n".join(lines)
REPORT.write_text(text + "\n")
print("\n" + text)
print(f"\nWrote {REPORT}")
