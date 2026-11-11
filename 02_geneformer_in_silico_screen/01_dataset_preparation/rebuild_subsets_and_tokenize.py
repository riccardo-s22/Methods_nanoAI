#!/usr/bin/env python
"""Rebuild cell-type subsets from the paper-QC AnnData and refresh the
geneformer_inputs/ folders. Tokenization is done by a separate script that
runs in the `geneformer` env (this one runs in `geneformer_prep`)."""
import os
import shutil
from pathlib import Path
import numpy as np
import scanpy as sc
from scipy import sparse

# Section path configuration (GF_WORK, GF_GEO, ...): see ../gf_paths.py.
import os, sys  # noqa: E401
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gf_paths  # noqa: E402

BASE = Path(gf_paths.GF_WORK)
IN = BASE / "h5ad/Garza_GSE209552_all_raw_counts_paperQC_annotated.h5ad"
OUTDIR = BASE / "h5ad/celltype_subsets"
INPUTS = BASE / "h5ad/geneformer_inputs"

adata = sc.read_h5ad(IN)
for col in ["ensembl_id"]:
    assert col in adata.var.columns
for col in ["n_counts", "condition", "cell_type", "filter_pass"]:
    assert col in adata.obs.columns

if not sparse.issparse(adata.X):
    adata.X = sparse.csr_matrix(adata.X)

adata.obs["condition"] = adata.obs["condition"].astype(str).replace({
    "ctrl": "control", "Control": "control", "tbi": "TBI",
})

SUBSETS = {
    "oligodendrocytes_raw_geneformer_ready.h5ad":     ["oligodendrocyte"],
    "OPCs_raw_geneformer_ready.h5ad":                 ["OPC"],
    "excitatory_neurons_raw_geneformer_ready.h5ad":   ["excitatory_neuron"],
    "inhibitory_neurons_raw_geneformer_ready.h5ad":   ["inhibitory_neuron"],
    "astrocytes_raw_geneformer_ready.h5ad":           ["astrocyte"],
    "microglia_macrophages_raw_geneformer_ready.h5ad":["microglia_macrophage"],
}
INPUT_MAP = {
    "oligo":               "oligodendrocytes_raw_geneformer_ready.h5ad",
    "OPC":                 "OPCs_raw_geneformer_ready.h5ad",
    "excitatory_neuron":   "excitatory_neurons_raw_geneformer_ready.h5ad",
    "inhibitory_neuron":   "inhibitory_neurons_raw_geneformer_ready.h5ad",
    "astrocyte":           "astrocytes_raw_geneformer_ready.h5ad",
    "microglia_macrophage":"microglia_macrophages_raw_geneformer_ready.h5ad",
}

# Wipe stale outputs.
for d in (OUTDIR, INPUTS):
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True, exist_ok=True)

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
    lines.append(f"WROTE {fname}: {sub.n_obs} cells x {sub.n_vars} genes  conditions={ct}  "
                 f"n_samples={sub.obs['sample'].nunique()}")

# Populate per-subset Geneformer input folders (hardlink for speed on NTFS).
for sub_dir, fname in INPUT_MAP.items():
    d = INPUTS / sub_dir
    d.mkdir(parents=True, exist_ok=True)
    src = OUTDIR / fname
    dst = d / fname
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)
    lines.append(f"  linked geneformer_inputs/{sub_dir}/{fname}")

(OUTDIR / "subset_creation_report.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
