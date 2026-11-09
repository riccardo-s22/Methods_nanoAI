#!/usr/bin/env python
"""Validate Geneformer-ready .h5ad files against the tokenizer contract."""
from pathlib import Path
import numpy as np
import scanpy as sc
from scipy import sparse

BASE = Path("/mnt/d/TBI_SC/Geneformer_TBI")
FILES = sorted(list((BASE / "h5ad").glob("*.h5ad")) + list((BASE / "h5ad/celltype_subsets").glob("*.h5ad")))

for f in FILES:
    print("\n" + "=" * 100)
    print(f)
    ad = sc.read_h5ad(f)
    print(ad)
    print("X sparse:", sparse.issparse(ad.X), "| X dtype:", ad.X.dtype)

    checks = {
        "var.ensembl_id":     "ensembl_id" in ad.var.columns,
        "obs.n_counts":       "n_counts" in ad.obs.columns,
        "obs.filter_pass":    "filter_pass" in ad.obs.columns,
        "obs.condition":      "condition" in ad.obs.columns,
        "obs.cell_type":      "cell_type" in ad.obs.columns,
    }
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")

    if "ensembl_id" in ad.var.columns:
        ens = ad.var["ensembl_id"].astype(str)
        frac = float(ens.str.startswith("ENSG").mean())
        print(f"  ensembl_id ENSG-like fraction: {frac:.3f}  examples: {ens.head(3).tolist()}")

    if "n_counts" in ad.obs.columns:
        s = ad.obs["n_counts"].astype(float)
        print(f"  n_counts: min={s.min():.0f}  median={s.median():.0f}  max={s.max():.0f}")

    if "condition" in ad.obs.columns:
        print("  condition counts:", ad.obs["condition"].value_counts(dropna=False).to_dict())

    if "cell_type" in ad.obs.columns:
        print("  cell_type counts:", ad.obs["cell_type"].value_counts(dropna=False).to_dict())

    x = ad.X.data if sparse.issparse(ad.X) else np.asarray(ad.X).ravel()
    x = x[:min(len(x), 200000)]
    if len(x):
        print(f"  X data: min={x.min()} max={x.max()} integer-like={bool(np.allclose(x, np.round(x)))}")

    print(f"  unique cell IDs: {ad.obs_names.is_unique}")

print("\nValidation complete.")
