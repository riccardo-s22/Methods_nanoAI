#!/usr/bin/env python
"""Apply the paper's additional hard cutoff (n_genes >= 1000) on top of the
paper-style MAD QC + scrublet filter, save a paper-matched primary AnnData,
and write a focused set of cell-type subsets.

Reports:
  - total cells, TBI/control
  - cells per sample, per cell type, per (cell_type x condition)
  - genes detected in excitatory neurons (>=1 read in >=1 cell)
  - overlap of those genes with the Geneformer V2 token dictionary
    (token_dictionary_gc104M.pkl)
"""
import pickle
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

BASE = Path("/mnt/d/TBI_SC/Geneformer_TBI")
IN = BASE / "h5ad/Garza_GSE209552_all_raw_counts_paperQC_annotated.h5ad"
OUT = BASE / "h5ad/Garza_GSE209552_all_raw_counts_paperQC_nGene1000_annotated.h5ad"
SUBSET_DIR = BASE / "h5ad/celltype_subsets_paperQC_nGene1000"
SUBSET_DIR.mkdir(parents=True, exist_ok=True)
REPORT = BASE / "reports/paperQC_nGene1000_report.txt"
V2_TOKEN_PKL = Path("/home/riccardo-s/miniconda3/envs/geneformer/lib/python3.10/site-packages/geneformer/token_dictionary_gc104M.pkl")

PAPER_TOTAL = 15110
PAPER_TBI = 6806
PAPER_CONTROL = 8304

# ---------------------------------------------------------------------
print(f"Reading {IN.name} ...", flush=True)
adata = sc.read_h5ad(IN)
print(f"  Input (MAD QC + doublet-filtered): {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

# Recompute n_genes_by_counts from current .X to be 100% explicit
# (.X is raw counts; n_genes_by_counts may have been computed on the pre-QC matrix).
n_genes_recalc = np.asarray((adata.X > 0).sum(axis=1)).ravel()
adata.obs["n_genes_by_counts_recalc"] = n_genes_recalc

# Hard floor.
mask = n_genes_recalc >= 1000
print(f"  cells with n_genes >= 1000: {int(mask.sum())} / {adata.n_obs}", flush=True)

final = adata[mask].copy()
final.obs["n_counts"] = np.asarray(final.X.sum(axis=1)).ravel().astype(np.float32)
final.obs["filter_pass"] = 1

# Drop the "removed_qc" cell_type bucket (would be empty since input is QC-pass,
# but be defensive).
final = final[final.obs["cell_type"].astype(str) != "removed_qc"].copy()

assert sparse.issparse(final.X)
x = final.X.data[:200000]
assert np.allclose(x, np.round(x)), "Raw counts must be integer-like."

final.write_h5ad(OUT, compression="gzip")
print(f"\nWrote {OUT}  ({final.n_obs} cells x {final.n_vars} genes)", flush=True)

# ---------------------------------------------------------------------
# Cell-type subsets requested.
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
    out_path = SUBSET_DIR / fname
    sub.write_h5ad(out_path, compression="gzip")
    cc = sub.obs["condition"].value_counts(dropna=False).to_dict()
    subset_info[fname] = (sub.n_obs, cc)
    print(f"  wrote {fname}: {sub.n_obs} cells  conditions={cc}", flush=True)

# ---------------------------------------------------------------------
# Genes detected in excitatory neurons + V2 vocab overlap.
exc = final[final.obs["cell_type"].astype(str) == "excitatory_neuron"]
exc_X = exc.X
gene_present = np.asarray((exc_X > 0).sum(axis=0)).ravel()
detected_mask = gene_present >= 1
exc_detected_ens = final.var["ensembl_id"].astype(str).values[detected_mask]
n_exc_detected = int(detected_mask.sum())

with open(V2_TOKEN_PKL, "rb") as fh:
    token_dict = pickle.load(fh)
v2_ens = {k for k in token_dict if isinstance(k, str) and k.startswith("ENSG")}
in_vocab = set(exc_detected_ens) & v2_ens
n_in_vocab = len(in_vocab)

# ---------------------------------------------------------------------
# Report.
lines = []
lines.append("=" * 70)
lines.append("Paper-matched primary object: paperQC + scrublet + n_genes>=1000")
lines.append("=" * 70)
lines.append(f"Total cells:           {final.n_obs}  (paper: {PAPER_TOTAL})")
tbi_n = int((final.obs['condition'] == 'TBI').sum())
ctrl_n = int((final.obs['condition'] == 'control').sum())
lines.append(f"TBI cells:             {tbi_n}  (paper: {PAPER_TBI})")
lines.append(f"control cells:         {ctrl_n}  (paper: {PAPER_CONTROL})")
lines.append("")
lines.append("Cells per sample:")
lines.append(final.obs["sample"].value_counts().sort_index().to_string())
lines.append("")
lines.append("Cells per cell_type:")
ct_counts = final.obs["cell_type"].value_counts(dropna=False)
lines.append(ct_counts.to_string())
lines.append("")
lines.append("Cells per (cell_type x condition):")
lines.append(pd.crosstab(final.obs["cell_type"], final.obs["condition"]).to_string())
lines.append("")
lines.append("Core cell-type counts requested:")
for ct in ["excitatory_neuron", "inhibitory_neuron", "oligodendrocyte", "OPC", "astrocyte", "microglia_macrophage"]:
    n = int((final.obs["cell_type"] == ct).sum())
    lines.append(f"  {ct:24s} {n}")
lines.append("")
lines.append("Excitatory neuron gene counts:")
lines.append(f"  Total genes in matrix:               {final.n_vars}")
lines.append(f"  Genes detected (>=1 read >=1 cell):  {n_exc_detected}")
lines.append(f"  Geneformer V2 vocab size (ENSG):     {len(v2_ens)}")
lines.append(f"  Excitatory-detected in V2 vocab:     {n_in_vocab}")
lines.append("")
lines.append("Subset h5ad outputs:")
for fname, (n, cc) in subset_info.items():
    lines.append(f"  {fname}: n={n}, conditions={cc}")
lines.append("")
lines.append(f"Primary object: {OUT}")
report_text = "\n".join(lines)
REPORT.write_text(report_text + "\n")
print("\n" + report_text)
print(f"\nWrote {REPORT}")
