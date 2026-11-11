#!/usr/bin/env python
"""Carve cell-type-specific Geneformer-ready .h5ad files from the annotated AnnData."""
from pathlib import Path
import numpy as np
import scanpy as sc
from scipy import sparse

# Section path configuration (GF_WORK, GF_GEO, ...): see ../gf_paths.py.
import os, sys  # noqa: E401
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gf_paths  # noqa: E402

BASE = Path(gf_paths.GF_WORK)
IN = BASE / "h5ad/Garza_GSE209552_all_raw_counts_annotated_geneformer_ready.h5ad"
OUTDIR = BASE / "h5ad/celltype_subsets"
OUTDIR.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(IN)

for col in ["ensembl_id"]:
    assert col in adata.var.columns, f"Missing adata.var['{col}']"
for col in ["n_counts", "condition", "cell_type"]:
    assert col in adata.obs.columns, f"Missing adata.obs['{col}']"

if not sparse.issparse(adata.X):
    adata.X = sparse.csr_matrix(adata.X)

adata.obs["condition"] = adata.obs["condition"].astype(str).replace({
    "ctrl": "control", "Control": "control", "CONTROL": "control",
    "tbi": "TBI", "Traumatic brain injury": "TBI",
})

# Restrict to QC-passing cells so each subset is clean for Geneformer.
qc_pass = adata.obs.get("passed_qc", 1)
if hasattr(qc_pass, "astype"):
    adata = adata[adata.obs["passed_qc"].astype(int) == 1].copy()

SUBSETS = {
    "oligodendrocytes_raw_geneformer_ready.h5ad": ["oligodendrocyte"],
    "OPCs_raw_geneformer_ready.h5ad":             ["OPC"],
    "excitatory_neurons_raw_geneformer_ready.h5ad": ["excitatory_neuron"],
    "inhibitory_neurons_raw_geneformer_ready.h5ad": ["inhibitory_neuron"],
    "astrocytes_raw_geneformer_ready.h5ad":         ["astrocyte"],
    "microglia_macrophages_raw_geneformer_ready.h5ad": ["microglia_macrophage"],
}

lines = []
for fname, labels in SUBSETS.items():
    mask = adata.obs["cell_type"].astype(str).isin(labels)
    sub = adata[mask].copy()
    if sub.n_obs == 0:
        lines.append(f"SKIP {fname}: 0 cells for {labels}")
        continue
    sub.obs["n_counts"] = np.asarray(sub.X.sum(axis=1)).ravel().astype(np.float32)
    sub.obs["filter_pass"] = 1
    out = OUTDIR / fname
    sub.write_h5ad(out, compression="gzip")
    ct = sub.obs["condition"].value_counts(dropna=False).to_dict()
    n_samples = sub.obs["sample"].nunique() if "sample" in sub.obs.columns else "NA"
    lines.append(f"WROTE {fname}: {sub.n_obs} cells x {sub.n_vars} genes  conditions={ct}  n_samples={n_samples}")

report = OUTDIR / "subset_creation_report.txt"
report.write_text("\n".join(lines) + "\n")
print(report.read_text())
