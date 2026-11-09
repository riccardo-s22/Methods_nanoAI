#!/usr/bin/env python
"""De novo cell-type annotation for the merged GSE209552 snRNA-seq AnnData.

Strategy
--------
1. Read the raw-counts merged h5ad (raw UMI in .X).
2. Keep raw counts intact; do all clustering on a normalized COPY.
3. Apply light QC: n_counts >= 500, pct_mt < 20.
4. Normalize + log1p, HVG, PCA, kNN, Leiden.
5. Score canonical brain cell-type marker sets with sc.tl.score_genes.
6. Per cluster, pick the highest-scoring cell type; mark "unspecified" if no clear winner.
7. Write cell_type back to the FULL (un-QC'd) AnnData so we do not lose cells from
   the Geneformer-ready file; cells removed by QC get cell_type = "unspecified".
8. Save Garza_GSE209552_all_raw_counts_annotated_geneformer_ready.h5ad with raw counts in .X.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

BASE = Path("/mnt/d/TBI_SC/Geneformer_TBI")
IN = BASE / "h5ad/Garza_GSE209552_all_raw_counts_geneformer_ready.h5ad"
OUT = BASE / "h5ad/Garza_GSE209552_all_raw_counts_annotated_geneformer_ready.h5ad"
REPORT = BASE / "reports/cell_type_annotation_report.txt"

MARKERS = {
    "oligodendrocyte":      ["PLP1", "MBP", "MOG", "MAG", "MOBP", "CLDN11", "OPALIN"],
    "OPC":                  ["VCAN", "PDGFRA", "CSPG4", "SOX10", "OLIG1", "OLIG2"],
    "excitatory_neuron":    ["RBFOX3", "SLC17A7", "SATB2", "CAMK2A", "NRGN", "NEUROD2", "NEUROD6"],
    "inhibitory_neuron":    ["GAD1", "GAD2", "SLC6A1", "SLC32A1", "PVALB", "SST", "VIP"],
    "astrocyte":            ["GFAP", "AQP4", "SLC1A2", "ALDH1L1", "SLC1A3", "S100B"],
    "microglia_macrophage": ["FYB1", "PTPRC", "C3", "CX3CR1", "AIF1", "CSF1R", "P2RY12", "TMEM119"],
    "endothelial":          ["CLDN5", "FLT1", "PECAM1", "VWF", "CDH5"],
}

sc.settings.verbosity = 2

print(f"Reading {IN} ...", flush=True)
adata_full = sc.read_h5ad(IN)
print(adata_full, flush=True)

# Confirm raw counts are in .X.
assert sparse.issparse(adata_full.X) or hasattr(adata_full.X, "shape")
print("X dtype:", adata_full.X.dtype, "| sparse:", sparse.issparse(adata_full.X), flush=True)

# Build a gene_symbol -> ensembl_id reverse map for marker lookup.
if "gene_symbol" not in adata_full.var.columns:
    raise SystemExit("Need gene_symbol in adata.var to look up markers.")
sym2ens = {s: e for e, s in zip(adata_full.var["ensembl_id"].astype(str), adata_full.var["gene_symbol"].astype(str)) if s}

def markers_to_ens(symbols):
    return [sym2ens[s] for s in symbols if s in sym2ens]

marker_ens = {ct: markers_to_ens(syms) for ct, syms in MARKERS.items()}
print("Marker counts:", {k: len(v) for k, v in marker_ens.items()}, flush=True)

# ---- QC + clustering on a copy ----
qc = adata_full.copy()
qc.var["mt"] = qc.var["gene_symbol"].astype(str).str.startswith("MT-")
sc.pp.calculate_qc_metrics(qc, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)

# Light filter to remove debris/dying nuclei before clustering.
pre = qc.n_obs
qc = qc[(qc.obs["total_counts"] >= 500) & (qc.obs["pct_counts_mt"] < 20)].copy()
print(f"QC-pass cells: {qc.n_obs} / {pre}", flush=True)

# Normalize + log on the QC copy. Keep adata_full.X untouched.
sc.pp.normalize_total(qc, target_sum=1e4)
sc.pp.log1p(qc)

sc.pp.highly_variable_genes(qc, n_top_genes=3000, flavor="seurat_v3" if False else "seurat", batch_key=None)
sc.pp.scale(qc, max_value=10)
sc.tl.pca(qc, n_comps=50)
sc.pp.neighbors(qc, n_neighbors=15, n_pcs=30)
sc.tl.leiden(qc, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
print("Cluster counts:")
print(qc.obs["leiden"].value_counts().sort_index().to_string(), flush=True)

# Score canonical markers per cell using log-normalized data.
for ct, ens_list in marker_ens.items():
    if not ens_list:
        qc.obs[f"score_{ct}"] = 0.0
        continue
    sc.tl.score_genes(qc, gene_list=ens_list, score_name=f"score_{ct}", use_raw=False)

# Mean score per cluster → assign cell type.
score_cols = [f"score_{ct}" for ct in MARKERS]
cluster_means = qc.obs.groupby("leiden", observed=False)[score_cols].mean()
cluster_means.columns = list(MARKERS.keys())

# Assign cluster -> cell type by argmax, but require a minimum score margin to avoid spurious calls.
def assign_row(row, margin=0.05):
    top = row.sort_values(ascending=False)
    if top.iloc[0] - top.iloc[1] < margin:
        return "unspecified"
    return top.index[0]

cluster_to_ct = cluster_means.apply(assign_row, axis=1)
print("\nCluster -> cell_type:")
print(pd.concat([cluster_to_ct.rename("cell_type"), cluster_means.round(3)], axis=1).to_string(), flush=True)

qc.obs["cell_type"] = qc.obs["leiden"].map(cluster_to_ct).astype(str)

# Propagate cell_type back to the FULL adata (cells removed by QC -> 'unspecified').
ct_full = pd.Series("unspecified", index=adata_full.obs_names)
ct_full.loc[qc.obs_names] = qc.obs["cell_type"].values
adata_full.obs["cell_type"] = ct_full.values
adata_full.obs["passed_qc"] = adata_full.obs_names.isin(qc.obs_names).astype(int)
# Also persist leiden and per-cell scores for traceability.
adata_full.obs["leiden"] = pd.Series(qc.obs["leiden"]).reindex(adata_full.obs_names).astype("object").fillna("NA").values
for c in score_cols:
    adata_full.obs[c] = pd.Series(qc.obs[c]).reindex(adata_full.obs_names).astype(float).fillna(np.nan).values

# Sanity: raw counts preserved.
assert sparse.issparse(adata_full.X)
x = adata_full.X.data[:100000]
assert np.allclose(x, np.round(x)), "Raw counts in .X should remain integer-like."

adata_full.write_h5ad(OUT, compression="gzip")
print(f"\nWrote {OUT}", flush=True)

with open(REPORT, "w") as fh:
    fh.write(str(adata_full) + "\n\n")
    fh.write("QC filter: total_counts >= 500 and pct_counts_mt < 20\n")
    fh.write(f"Cells kept by QC: {qc.n_obs} / {pre}\n\n")
    fh.write("Cluster -> cell_type assignment (rows=Leiden cluster):\n")
    fh.write(pd.concat([cluster_to_ct.rename("cell_type"), cluster_means.round(3)], axis=1).to_string() + "\n\n")
    fh.write("Cell type counts (full data, including QC-failed -> 'unspecified'):\n")
    fh.write(adata_full.obs["cell_type"].value_counts(dropna=False).to_string() + "\n\n")
    fh.write("Cell type x condition (full data):\n")
    fh.write(pd.crosstab(adata_full.obs["cell_type"], adata_full.obs["condition"]).to_string() + "\n")
print(f"Wrote {REPORT}")
