#!/usr/bin/env python
"""Gene-expression UMAP of ALL cell types in the paperQC + nGene1000 primary
object (15,172 cells). Standard scanpy workflow on raw counts:

    normalize_total(1e4) -> log1p -> 2,000 HVGs (seurat, batch=sample)
    -> scale(max 10) -> PCA(50) -> kNN(30 neighbors, 30 PCs) -> UMAP

No batch correction (raw donor/sample structure is left visible on purpose so
donor effects can be judged in the donor panel).

Panels: cell_type, condition, donor, region, n_genes (log10), pct_counts_mt.

Output: results/umap/all_celltypes_expr_umap.png (+ .pdf, coords .tsv.gz)
"""
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[1]
H5AD = BASE / "h5ad/Garza_GSE209552_all_raw_counts_paperQC_nGene1000_annotated.h5ad"
OUT_DIR = BASE / "results/umap"
PREFIX = "all_celltypes_expr"
SEED = 0

CT_ORDER = ["excitatory_neuron", "inhibitory_neuron", "oligodendrocyte", "OPC",
            "astrocyte", "microglia_macrophage", "endothelial",
            "pericyte_mural", "unspecified"]
CT_COLORS = {
    "excitatory_neuron": "#1f77b4", "inhibitory_neuron": "#17becf",
    "oligodendrocyte": "#2ca02c", "OPC": "#bcbd22", "astrocyte": "#ff7f0e",
    "microglia_macrophage": "#d62728", "endothelial": "#9467bd",
    "pericyte_mural": "#e377c2", "unspecified": "#7f7f7f",
}


def scatter_cat(ax, df, col, order=None, colors=None, title=None):
    cats = order or sorted(df[col].astype(str).unique(),
                           key=lambda s: (len(s), s))
    cmap = plt.cm.tab20
    for i, c in enumerate(cats):
        m = df[col].astype(str) == c
        if not m.any():
            continue
        color = colors[c] if colors else cmap(i % 20)
        ax.scatter(df.loc[m, "UMAP1"], df.loc[m, "UMAP2"], s=2, alpha=0.6,
                   linewidths=0, color=color, rasterized=True,
                   label=f"{c} (n={int(m.sum()):,})")
    ax.set_title(title or f"by {col}")
    ax.legend(markerscale=5, fontsize=7, loc="center left",
              bbox_to_anchor=(1.0, 0.5), frameon=False)


def scatter_cont(ax, fig, df, col, title):
    o = np.argsort(df[col].to_numpy())
    sca = ax.scatter(df["UMAP1"].iloc[o], df["UMAP2"].iloc[o], c=df[col].iloc[o],
                     s=2, linewidths=0, cmap="viridis", rasterized=True)
    fig.colorbar(sca, ax=ax, shrink=0.8)
    ax.set_title(title)


def label_centroids(ax, df):
    for ct, g in df.groupby("cell_type", observed=True):
        if len(g) < 50:
            continue
        ax.text(g["UMAP1"].median(), g["UMAP2"].median(), ct, fontsize=8,
                ha="center", va="center", weight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7,
                          lw=0))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(H5AD)
    print(adata)
    adata.var_names = adata.var["gene_symbol"].astype(str)
    adata.var_names_make_unique()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat",
                                batch_key="sample")
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, random_state=SEED)
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30, random_state=SEED)
    sc.tl.umap(adata, min_dist=0.3, random_state=SEED)

    df = adata.obs[["cell_type", "condition", "donor", "region", "sample",
                    "n_genes", "n_counts", "pct_counts_mt"]].copy()
    df["UMAP1"], df["UMAP2"] = adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1]
    df["log10_n_genes"] = np.log10(df["n_genes"].astype(float))
    df.to_csv(OUT_DIR / f"{PREFIX}_umap_coords.tsv.gz", sep="\t")

    fig, axes = plt.subplots(2, 3, figsize=(21, 12))
    ax = axes[0, 0]
    scatter_cat(ax, df, "cell_type", order=CT_ORDER, colors=CT_COLORS,
                title="Cell type")
    label_centroids(ax, df)
    scatter_cat(axes[0, 1], df, "condition", order=["control", "TBI"],
                colors={"control": "#4575b4", "TBI": "#d73027"},
                title="Condition")
    scatter_cat(axes[0, 2], df, "donor", title="Donor")
    scatter_cat(axes[1, 0], df, "region", title="Region")
    scatter_cont(axes[1, 1], fig, df, "log10_n_genes", "log10(n_genes)")
    scatter_cont(axes[1, 2], fig, df, "pct_counts_mt", "% mitochondrial counts")
    for a in axes.flat:
        a.set_xticks([]); a.set_yticks([])
        a.set_xlabel("UMAP1"); a.set_ylabel("UMAP2")
    fig.suptitle("Garza GSE209552 — all cell types, gene-expression UMAP "
                 f"(paperQC + nGene≥1000, n={len(df):,} cells; "
                 "2,000 HVGs, 30 PCs, no batch correction)", fontsize=13)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = OUT_DIR / f"{PREFIX}_umap.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
