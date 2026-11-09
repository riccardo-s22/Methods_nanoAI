#!/usr/bin/env python
"""Rebuild the annotated AnnData using Garza et al. 2023's QC pipeline.

Issues addressed (vs the initial build):
  (1) Lax fixed-threshold QC -> paper-style per-sample MAD filter (3 MADs) on
      pct_mt and n_genes, matching src/sc_pipeline/do_quality_filtering.R.
  (4) Cell-type labels applied to a too-permissive matrix -> annotate ONLY on
      cells that pass the paper's filter + a scrublet doublet flag, on a
      Harmony-integrated embedding (paper uses Harmony).

Inputs:  h5ad/Garza_GSE209552_all_raw_counts_geneformer_ready.h5ad (raw merged)
Outputs:
  - h5ad/Garza_GSE209552_all_raw_counts_paperQC_annotated.h5ad
      Annotated AnnData restricted to QC-passing, doublet-negative cells.
      .X = raw UMI counts. Adds: cell_type, leiden_paperQC, doublet_score,
                                  predicted_doublet, qc_pass_mad.
  - reports/paperQC_rebuild_report.txt
"""
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from scipy.stats import median_abs_deviation
import harmonypy as hm
import scrublet as scr

BASE = Path("/mnt/d/TBI_SC/Geneformer_TBI")
IN = BASE / "h5ad/Garza_GSE209552_all_raw_counts_geneformer_ready.h5ad"
OUT = BASE / "h5ad/Garza_GSE209552_all_raw_counts_paperQC_annotated.h5ad"
REPORT = BASE / "reports/paperQC_rebuild_report.txt"

# Canonical brain marker panels (gene symbols). Includes paper-cluster-derived
# markers from data/MarkersPerCluster.tsv (Garza 2023).
MARKERS = {
    "oligodendrocyte":    ["PLP1", "MBP", "MOG", "MAG", "MOBP", "CLDN11", "OPALIN",
                           "ST18", "RNF220", "SLC44A1", "CTNNA3"],
    "OPC":                ["VCAN", "PDGFRA", "CSPG4", "SOX10", "OLIG1", "OLIG2",
                           "LHFPL3", "PTPRZ1", "PCDH15", "TNR", "MEGF11"],
    "excitatory_neuron":  ["SLC17A7", "SATB2", "CAMK2A", "NRGN", "NEUROD2", "NEUROD6",
                           "CBLN2", "CUX2", "HOMER1", "FAM19A1", "TSHZ2", "FOXP2"],
    "inhibitory_neuron":  ["GAD1", "GAD2", "SLC6A1", "SLC32A1", "PVALB", "SST", "VIP",
                           "ADARB2", "DLX6-AS1", "NXPH1", "SYNPR"],
    "astrocyte":          ["GFAP", "AQP4", "SLC1A2", "ALDH1L1", "SLC1A3", "S100B",
                           "ADGRV1", "GPC5"],
    "microglia_macrophage":["PTPRC", "C3", "CX3CR1", "AIF1", "CSF1R", "P2RY12", "TMEM119",
                            "SPP1", "APBB1IP", "DOCK8", "RUNX1", "FYB1"],
    "endothelial":        ["CLDN5", "FLT1", "PECAM1", "VWF", "CDH5", "MECOM", "ABCB1"],
    "pericyte_mural":     ["LAMA2", "EBF1", "DLC1", "PDGFRB", "RGS5", "ANPEP"],
}

sc.settings.verbosity = 1

print(f"Reading {IN.name} ...", flush=True)
adata = sc.read_h5ad(IN)
print(adata, flush=True)
n0 = adata.n_obs

# Build symbol -> ensembl_id map for marker scoring.
sym2ens = {s: e for e, s in zip(adata.var["ensembl_id"].astype(str),
                                adata.var["gene_symbol"].astype(str)) if s}
def markers_to_ens(syms):
    return [sym2ens[s] for s in syms if s in sym2ens]
marker_ens = {ct: markers_to_ens(syms) for ct, syms in MARKERS.items()}

# Ensure QC metrics exist.
if "pct_counts_mt" not in adata.obs.columns:
    adata.var["mt"] = adata.var["gene_symbol"].astype(str).str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)

# -----------------------------------------------------------------
# Paper-style QC: per-sample MAD filter (3 MAD pct_mt, 3 MAD n_genes lower tail).
# Matches src/sc_pipeline/do_quality_filtering.R defaults FILTER_MITO_NMADS=3,
# FILTER_LOW_NUM_GENES_NMADS=3.
# -----------------------------------------------------------------
print("\nApplying paper-style per-sample MAD QC (3 MAD pct_mt high, 3 MAD n_genes low) ...", flush=True)
NMADS = 3.0
keep = np.zeros(adata.n_obs, dtype=bool)
qc_lines = []
for s, idx in adata.obs.groupby("sample", observed=False).groups.items():
    pct_mt = adata.obs.loc[idx, "pct_counts_mt"].astype(float).values
    n_genes = adata.obs.loc[idx, "n_genes_by_counts"].astype(float).values
    med_mt = np.median(pct_mt);  mad_mt = median_abs_deviation(pct_mt) or 1e-6
    med_g  = np.median(n_genes); mad_g  = median_abs_deviation(n_genes) or 1e-6
    pass_mt = pct_mt <= med_mt + NMADS * mad_mt
    pass_g  = n_genes >= med_g - NMADS * mad_g
    # Hard floor: drop empty/near-empty droplets (n_genes < 200).
    pass_g  = pass_g & (n_genes >= 200)
    pass_all = pass_mt & pass_g
    sample_keep = pd.Series(False, index=adata.obs_names)
    sample_keep.loc[idx] = pass_all
    keep |= sample_keep.values
    qc_lines.append(f"  {s}: n={len(idx)}  median_mt={med_mt:.2f}±{mad_mt:.2f}  "
                    f"median_genes={med_g:.0f}±{mad_g:.0f}  pass={int(pass_all.sum())}")
adata.obs["qc_pass_mad"] = keep.astype(int)
print("\n".join(qc_lines), flush=True)
print(f"  Total kept by MAD QC: {keep.sum()} / {n0}", flush=True)

# -----------------------------------------------------------------
# Doublet detection with Scrublet (per sample) on raw counts.
# -----------------------------------------------------------------
print("\nScrublet per sample ...", flush=True)
doublet_score = np.full(adata.n_obs, np.nan)
predicted_doublet = np.zeros(adata.n_obs, dtype=bool)
for s, idx in adata.obs.groupby("sample", observed=False).groups.items():
    sub_pos = np.where(adata.obs_names.isin(idx))[0]
    X = adata.X[sub_pos]
    try:
        s_scr = scr.Scrublet(X, expected_doublet_rate=0.08)
        scores, called = s_scr.scrub_doublets(min_counts=1, min_cells=3,
                                              min_gene_variability_pctl=85,
                                              n_prin_comps=30,
                                              verbose=False)
        if called is None:
            called = scores > 0.25
    except Exception as e:
        print(f"  scrublet failed on {s}: {e}; defaulting to no doublet flag")
        scores = np.zeros(len(sub_pos))
        called = np.zeros(len(sub_pos), dtype=bool)
    doublet_score[sub_pos] = scores
    predicted_doublet[sub_pos] = called
    print(f"  {s}: doublets called {int(called.sum())}/{len(sub_pos)}", flush=True)
adata.obs["doublet_score"] = doublet_score
adata.obs["predicted_doublet"] = predicted_doublet.astype(int)

# Combined pass = paper MAD QC AND not flagged as doublet.
adata.obs["qc_pass_final"] = ((adata.obs["qc_pass_mad"] == 1) &
                              (adata.obs["predicted_doublet"] == 0)).astype(int)
n_pass = int(adata.obs["qc_pass_final"].sum())
print(f"\nFinal QC pass (MAD + not-doublet): {n_pass} / {n0}", flush=True)

# -----------------------------------------------------------------
# Clustering on QC-pass cells with Harmony batch correction by sample.
# -----------------------------------------------------------------
print("\nClustering with Harmony batch correction ...", flush=True)
qc = adata[adata.obs["qc_pass_final"] == 1].copy()
qc.layers["counts"] = qc.X.copy()  # preserve raw counts for later
sc.pp.normalize_total(qc, target_sum=1e4)
sc.pp.log1p(qc)
sc.pp.highly_variable_genes(qc, n_top_genes=3000, flavor="seurat", batch_key="sample")
sc.pp.scale(qc, max_value=10)
sc.tl.pca(qc, n_comps=50)
# Harmony integrate by sample.
ho = hm.run_harmony(qc.obsm["X_pca"], qc.obs[["sample"]], "sample", max_iter_harmony=20)
Z = ho.Z_corr
qc.obsm["X_pca_harmony"] = Z.T if Z.shape[0] == qc.obsm["X_pca"].shape[1] else Z
sc.pp.neighbors(qc, n_neighbors=15, n_pcs=30, use_rep="X_pca_harmony")
sc.tl.leiden(qc, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
print("Clusters:")
print(qc.obs["leiden"].value_counts().sort_index().to_string(), flush=True)

# -----------------------------------------------------------------
# Per-cell marker scoring + per-cluster argmax for cell_type assignment.
# -----------------------------------------------------------------
for ct, ens_list in marker_ens.items():
    if ens_list:
        sc.tl.score_genes(qc, gene_list=ens_list, score_name=f"score_{ct}", use_raw=False)
    else:
        qc.obs[f"score_{ct}"] = 0.0
score_cols = [f"score_{ct}" for ct in MARKERS]
cluster_means = qc.obs.groupby("leiden", observed=False)[score_cols].mean()
cluster_means.columns = list(MARKERS.keys())

def assign_row(row, margin=0.05):
    top = row.sort_values(ascending=False)
    if top.iloc[0] - top.iloc[1] < margin:
        return "unspecified"
    return top.index[0]
cluster_to_ct = cluster_means.apply(assign_row, axis=1)
print("\nCluster -> cell_type (Harmony-integrated, paper-QC cells):")
print(pd.concat([cluster_to_ct.rename("cell_type"), cluster_means.round(3)], axis=1).to_string(), flush=True)
qc.obs["cell_type"] = qc.obs["leiden"].map(cluster_to_ct).astype(str)

# -----------------------------------------------------------------
# Propagate cell_type/leiden back to the FULL adata, then KEEP ONLY
# QC-pass cells for the output (concern #4: no labels on raw matrix).
# -----------------------------------------------------------------
ct_series = pd.Series("removed_qc", index=adata.obs_names)
ct_series.loc[qc.obs_names] = qc.obs["cell_type"].values
adata.obs["cell_type"] = ct_series.values
adata.obs["leiden_paperQC"] = pd.Series(qc.obs["leiden"]).reindex(adata.obs_names).astype("object").fillna("NA").values
for c in score_cols:
    adata.obs[c] = pd.Series(qc.obs[c]).reindex(adata.obs_names).astype(float).values

# Restrict outputs to QC-passing cells. This is the matrix the labels apply to.
final = adata[adata.obs["qc_pass_final"] == 1].copy()

# Sanity: raw integer counts intact.
assert sparse.issparse(final.X)
x = final.X.data[:200000]
assert np.allclose(x, np.round(x)), "Raw counts must be integer-like."

# Required Geneformer fields.
final.obs["n_counts"] = np.asarray(final.X.sum(axis=1)).ravel().astype(np.float32)
final.obs["filter_pass"] = 1

final.write_h5ad(OUT, compression="gzip")
print(f"\nWrote {OUT}  ({final.n_obs} cells x {final.n_vars} genes)", flush=True)

# Report.
lines = []
lines.append(f"Input cells: {n0}")
lines.append(f"Paper MAD QC pass: {int(adata.obs['qc_pass_mad'].sum())}")
lines.append(f"Doublets called: {int(adata.obs['predicted_doublet'].sum())}")
lines.append(f"Final (MAD + not-doublet): {n_pass}")
lines.append("")
lines.append("Per-sample QC stats:")
lines.extend(qc_lines)
lines.append("")
lines.append("Cluster -> cell_type assignment:")
lines.append(pd.concat([cluster_to_ct.rename("cell_type"), cluster_means.round(3)], axis=1).to_string())
lines.append("")
lines.append("Cell type counts (QC-passed cells):")
lines.append(final.obs["cell_type"].value_counts(dropna=False).to_string())
lines.append("")
lines.append("Cell type x condition:")
lines.append(pd.crosstab(final.obs["cell_type"], final.obs["condition"]).to_string())
lines.append("")
lines.append("Cells per sample (QC-passed):")
lines.append(final.obs["sample"].value_counts().sort_index().to_string())

REPORT.write_text("\n".join(lines))
print(f"Wrote {REPORT}")
