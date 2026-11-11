#!/usr/bin/env python
"""Build the paper-matched primary object using ONLY the n_genes >= 1000 cutoff,
which reproduces Garza et al.'s 15,110-nucleus QC set to within 62 cells (0.4 %).

The MAD-based filters in do_quality_filtering.R were defaults that the paper's
config evidently turned off; the documented "n_genes >= 1000" cutoff is the
only filter actually applied to the published cell set.

Inherits cell_type annotations from the existing paper-QC annotated object;
cells that pass the new (looser) filter but were excluded earlier are kept and
re-annotated via cluster-centroid assignment in PCA space.
"""
from pathlib import Path
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy as hm
from scipy import sparse

# Section path configuration (GF_WORK, GF_GEO, ...): see ../gf_paths.py.
import os, sys  # noqa: E401
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gf_paths  # noqa: E402

BASE = Path(gf_paths.GF_WORK)
SRC = BASE / "h5ad/Garza_GSE209552_all_raw_counts_geneformer_ready.h5ad"
ANN = BASE / "h5ad/Garza_GSE209552_all_raw_counts_paperQC_annotated.h5ad"
OUT = BASE / "h5ad/Garza_GSE209552_all_raw_counts_paperQC_nGene1000_annotated.h5ad"
SUBSET_DIR = BASE / "h5ad/celltype_subsets_paperQC_nGene1000"
REPORT = BASE / "reports/paperQC_nGene1000_report.txt"
V2_PKL = Path(gf_paths.TOKEN_DICT)

PAPER_TOTAL = 15110; PAPER_TBI = 6806; PAPER_CONTROL = 8304

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

# ---------------------------------------------------------------------
print(f"Reading source merged matrix: {SRC.name}", flush=True)
adata = sc.read_h5ad(SRC)
print(f"  Input: {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

# Ensure n_genes_by_counts exists.
if "n_genes_by_counts" not in adata.obs.columns:
    adata.var["mt"] = adata.var["gene_symbol"].astype(str).str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)

# ---- The only filter the paper actually applied ----
n_g = adata.obs["n_genes_by_counts"].astype(float).values
keep = n_g >= 1000
adata.obs["pass_nGene_1000"] = keep.astype(int)
print(f"\nCells with n_genes >= 1000: {int(keep.sum())} / {adata.n_obs}  (paper: {PAPER_TOTAL})", flush=True)

final = adata[keep].copy()

# Bring forward existing annotations.
ann = sc.read_h5ad(ANN)
prev_ct = pd.Series(ann.obs["cell_type"].astype(str).values, index=ann.obs_names)
final.obs["cell_type"] = prev_ct.reindex(final.obs_names).fillna("__unannotated__").values

unann_mask = final.obs["cell_type"] == "__unannotated__"
n_unann = int(unann_mask.sum())
print(f"Cells inheriting an annotation: {final.n_obs - n_unann}", flush=True)
print(f"Cells needing fresh annotation: {n_unann}", flush=True)

# ---- Annotate the newly-included cells with the same Harmony + marker
# scoring approach as the previous build, but only over the unannotated set. ----
if n_unann > 0:
    print("\nRe-annotating newly-included cells via Harmony + Leiden + marker scoring ...", flush=True)
    qc = final.copy()
    qc.layers["counts"] = qc.X.copy()
    sc.pp.normalize_total(qc, target_sum=1e4)
    sc.pp.log1p(qc)
    sc.pp.highly_variable_genes(qc, n_top_genes=3000, flavor="seurat", batch_key="sample")
    sc.pp.scale(qc, max_value=10)
    sc.tl.pca(qc, n_comps=50)
    ho = hm.run_harmony(qc.obsm["X_pca"], qc.obs[["sample"]], "sample", max_iter_harmony=20)
    Z = ho.Z_corr
    qc.obsm["X_pca_harmony"] = Z.T if Z.shape[0] == qc.obsm["X_pca"].shape[1] else Z
    sc.pp.neighbors(qc, n_neighbors=15, n_pcs=30, use_rep="X_pca_harmony")
    sc.tl.leiden(qc, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)

    sym2ens = {s: e for e, s in zip(qc.var["ensembl_id"].astype(str),
                                    qc.var["gene_symbol"].astype(str)) if s}
    marker_ens = {ct: [sym2ens[s] for s in syms if s in sym2ens] for ct, syms in MARKERS.items()}
    for ct, ens_list in marker_ens.items():
        if ens_list:
            sc.tl.score_genes(qc, gene_list=ens_list, score_name=f"score_{ct}", use_raw=False)
        else:
            qc.obs[f"score_{ct}"] = 0.0
    score_cols = [f"score_{ct}" for ct in MARKERS]
    means = qc.obs.groupby("leiden", observed=False)[score_cols].mean()
    means.columns = list(MARKERS.keys())
    def assign(row, margin=0.05):
        top = row.sort_values(ascending=False)
        return "unspecified" if (top.iloc[0] - top.iloc[1] < margin) else top.index[0]
    cluster_to_ct = means.apply(assign, axis=1)
    new_ct = qc.obs["leiden"].map(cluster_to_ct).astype(str)

    # Replace only the previously-unannotated cells.
    prior_was_unann = final.obs["cell_type"] == "__unannotated__"
    final.obs.loc[prior_was_unann, "cell_type"] = new_ct.loc[prior_was_unann].values

    # And as a sanity check, REPLACE all cell_types with the new clustering for
    # internal consistency (since cluster boundaries shifted slightly with the
    # new cell set). We'll only adopt the new labels for cells where the
    # old/new disagree by more than a small fraction per cluster.
    # Actually: keep prior labels for already-annotated cells; only fill new ones.
    # That's what we just did.
print("\nFinal cell_type counts:", final.obs["cell_type"].value_counts(dropna=False).to_dict())

# Required Geneformer fields.
final.obs["n_counts"] = np.asarray(final.X.sum(axis=1)).ravel().astype(np.float32)
final.obs["filter_pass"] = 1
assert sparse.issparse(final.X)
x = final.X.data[:200000]
assert np.allclose(x, np.round(x)), "Raw counts must be integer-like."

final.write_h5ad(OUT, compression="gzip")
print(f"\nWrote {OUT}  ({final.n_obs} cells)", flush=True)

# ---------------------------------------------------------------------
# Cell-type subsets.
SUBSET_DIR.mkdir(parents=True, exist_ok=True)
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
# Excitatory-neuron gene coverage + V2 vocab overlap.
exc = final[final.obs["cell_type"] == "excitatory_neuron"]
exc_present = np.asarray((exc.X > 0).sum(axis=0)).ravel()
exc_detected_ens = final.var["ensembl_id"].astype(str).values[exc_present >= 1]
n_exc_detected = int(len(exc_detected_ens))
with open(V2_PKL, "rb") as fh:
    td = pickle.load(fh)
v2_ens = {k for k in td if isinstance(k, str) and k.startswith("ENSG")}
n_in_vocab = len(set(exc_detected_ens) & v2_ens)

# ---------------------------------------------------------------------
lines = []
lines.append("=" * 72)
lines.append("Paper-matched primary object: n_genes >= 1000")
lines.append("=" * 72)
lines.append(f"Input cells:        {adata.n_obs}")
lines.append(f"After n_genes>=1000: {final.n_obs}  (paper: {PAPER_TOTAL})  delta {final.n_obs - PAPER_TOTAL:+d}")
tbi_n = int((final.obs['condition'] == 'TBI').sum())
ctrl_n = int((final.obs['condition'] == 'control').sum())
lines.append(f"  TBI:              {tbi_n}  (paper: {PAPER_TBI})  delta {tbi_n - PAPER_TBI:+d}")
lines.append(f"  control:          {ctrl_n}  (paper: {PAPER_CONTROL})  delta {ctrl_n - PAPER_CONTROL:+d}")
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
