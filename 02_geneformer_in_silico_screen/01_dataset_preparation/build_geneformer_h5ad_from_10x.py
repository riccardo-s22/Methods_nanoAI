#!/usr/bin/env python
"""Build a merged raw-counts .h5ad for Geneformer from the snRNA-seq GSE209552 samples.

Reads barcodes/features/matrix directly from /mnt/d/TBI_SC/GSE209552_RAW (NTFS — no symlinks)
to avoid the symlink-organization step on Windows-mounted drives.
"""
import gzip
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io as sio
from scipy import sparse

BASE = Path("/mnt/d/TBI_SC/Geneformer_TBI")
RAW = Path("/mnt/d/TBI_SC/GSE209552_RAW")
MANIFEST = BASE / "reports/snrna_sample_manifest.txt"
META = BASE / "metadata/gse209552_sample_metadata_curated.tsv"
OUTDIR = BASE / "h5ad"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT = OUTDIR / "Garza_GSE209552_all_raw_counts_geneformer_ready.h5ad"


def strip_ensembl_version(x: str) -> str:
    return re.sub(r"\.\d+$", "", str(x))


def read_one_sample(gsm: str, matrix_name: str, barcodes_name: str, features_name: str, sample_md: dict | None):
    matrix_path = RAW / matrix_name
    barcodes_path = RAW / barcodes_name
    features_path = RAW / features_name
    for p in (matrix_path, barcodes_path, features_path):
        if not p.exists():
            raise FileNotFoundError(p)

    # Matrix.mtx.gz: 10x writes genes x cells; transpose to cells x genes.
    with gzip.open(matrix_path, "rb") as fh:
        X = sio.mmread(fh).tocsr()  # genes x cells
    X = X.T.tocsr()  # cells x genes

    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t", names=["barcode"])
    # features.tsv.gz columns: ensembl_id, gene_symbol, feature_type
    features = pd.read_csv(features_path, header=None, sep="\t",
                           names=["ensembl_id", "gene_symbol", "feature_type"])
    # Geneformer is gene-expression only.
    gex_mask = features["feature_type"].astype(str).str.lower().str.startswith("gene")
    if not gex_mask.all():
        X = X[:, gex_mask.values]
        features = features.loc[gex_mask].reset_index(drop=True)

    features["ensembl_id"] = features["ensembl_id"].map(strip_ensembl_version)

    var = features.set_index("ensembl_id", drop=False)
    var.index.name = None
    # Make ensembl_id unique in case of duplicates
    if var.index.duplicated().any():
        dup_mask = var.index.duplicated(keep="first")
        keep = ~dup_mask
        X = X[:, keep.values]
        var = var.loc[keep]

    obs = pd.DataFrame({"barcode": barcodes["barcode"].values})
    obs["sample"] = gsm
    obs.index = [f"{gsm}_{bc}" for bc in obs["barcode"]]
    obs["cell_id"] = obs.index

    a = ad.AnnData(X=X.astype(np.float32), obs=obs, var=var)
    a.var_names_make_unique()

    # Inject sample-level metadata onto every cell.
    if sample_md is not None:
        for k, v in sample_md.items():
            if k in ("sample",):
                continue
            col = k if k not in a.obs.columns else f"samplemeta_{k}"
            a.obs[col] = v

    a.obs["condition"] = a.obs.get("condition", sample_md.get("condition", "unknown") if sample_md else "unknown")
    a.obs["assay"] = a.obs.get("assay", "snRNAseq_10x")
    a.obs["source"] = "GSE209552_GEO_processed_10x"
    a.obs["n_counts"] = np.asarray(a.X.sum(axis=1)).ravel()
    a.obs["n_genes"] = np.asarray((a.X > 0).sum(axis=1)).ravel()
    a.obs["filter_pass"] = 1
    return a


def main():
    if not MANIFEST.exists():
        sys.exit(f"Missing manifest: {MANIFEST}")
    manifest = [line.strip().split("|") for line in MANIFEST.read_text().splitlines() if line.strip()]
    md = pd.read_csv(META, sep="\t", dtype=str).fillna("").set_index("sample", drop=False)

    adatas = []
    for gsm, mtx, bc, ft in manifest:
        sm = md.loc[gsm].to_dict() if gsm in md.index else None
        a = read_one_sample(gsm, mtx, bc, ft, sm)
        cond = a.obs["condition"].iloc[0] if "condition" in a.obs.columns else "?"
        print(f"  {gsm}: {a.n_obs} cells x {a.n_vars} genes  condition={cond}", flush=True)
        adatas.append(a)

    print("Concatenating ...", flush=True)
    adata = ad.concat(adatas, axis=0, join="outer", merge="same", label="batch_sample", keys=[a.obs["sample"].iloc[0] for a in adatas], index_unique=None)

    # Re-attach .var ensembl_id properly (concat with outer join keeps the union of vars).
    # The var index after concat is the ensembl_id. Add explicit column.
    adata.var["ensembl_id"] = adata.var_names.astype(str)
    # Carry through gene_symbol from any subadata that has it.
    sym_map = {}
    for a in adatas:
        if "gene_symbol" in a.var.columns:
            sym_map.update(dict(zip(a.var.index.astype(str), a.var["gene_symbol"].astype(str))))
    if sym_map:
        adata.var["gene_symbol"] = pd.Series(adata.var_names.astype(str)).map(sym_map).fillna("").values

    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    adata.obs["filter_pass"] = 1
    adata.var["mt"] = adata.var.get("gene_symbol", pd.Series(adata.var_names, index=adata.var_names)).astype(str).str.startswith("MT-")

    try:
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
        if "total_counts" in adata.obs.columns:
            adata.obs["n_counts"] = adata.obs["total_counts"].astype(float)
    except Exception as e:
        print(f"QC warn: {e}")

    assert "ensembl_id" in adata.var.columns
    assert "n_counts" in adata.obs.columns
    assert sparse.issparse(adata.X)

    adata.write_h5ad(OUT, compression="gzip")
    print(f"Wrote {OUT}", flush=True)

    summary = []
    for col in ["condition", "sample", "assay", "region", "donor"]:
        if col in adata.obs.columns:
            summary.append((col, adata.obs[col].value_counts(dropna=False)))
    with open(OUTDIR / "h5ad_build_summary.txt", "w") as fh:
        fh.write(str(adata) + "\n\n")
        fh.write("obs cols:\n" + ", ".join(map(str, adata.obs.columns)) + "\n\n")
        fh.write("var cols:\n" + ", ".join(map(str, adata.var.columns)) + "\n\n")
        for col, s in summary:
            fh.write(f"-- {col} --\n{s.to_string()}\n\n")
    print("Done.")


if __name__ == "__main__":
    main()
