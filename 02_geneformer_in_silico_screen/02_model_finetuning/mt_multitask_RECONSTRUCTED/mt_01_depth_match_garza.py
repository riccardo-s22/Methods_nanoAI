"""Depth-match the Garza corpus for multi-task (MT) fine-tuning.

Quantile-bins the total counts per nucleus WITHIN each annotated cell type and
subsamples both conditions equally within each bin, so that the control-versus-TBI
sequencing-depth difference cannot act as a shortcut for the disease-state head.

Published effect of this step (Methods, "Multi-task fine-tuning"):
    15,034 -> 9,428 nuclei (4,714 control + 4,714 TBI)
    median control-vs-TBI depth difference reduced to <= 250 counts in every
    cell type, from 6,947 counts in excitatory neurons.

=============================================================================
Rewritten from the Methods after the original MT code was lost. It documents the
procedure but was not used to generate the published results.
=============================================================================
Values the Methods do not fix are marked UNRECORDED rather than guessed.

Original file name: 01b_depth_match_garza.py
=============================================================================
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# UNRECORDED: the number of quantile bins is not stated in the Methods and does
# not appear in the surviving notes. 10 deciles is the value most consistent with
# the reported outcome (<=250-count residual median gap at 9,428 retained nuclei);
# treat it as a reconstruction choice, not a recovered one.
N_BINS = 10

SEED = 42

CELL_TYPES = [
    "excitatory_neuron",
    "inhibitory_neuron",
    "oligodendrocyte",
    "OPC",
    "astrocyte",
    "microglia_macrophage",
]


def depth_match(obs: pd.DataFrame, n_bins: int = N_BINS, seed: int = SEED) -> np.ndarray:
    """Return the positional indices of the depth-matched nuclei.

    Within each cell type, nuclei are binned by quantiles of ``n_counts``
    (pooling both conditions so the bin edges are condition-independent). In
    each bin the two conditions are downsampled to the size of the smaller one,
    without replacement.
    """
    rng = np.random.default_rng(seed)
    keep: list[np.ndarray] = []

    for cell_type in sorted(obs["cell_type"].unique()):
        ct_mask = obs["cell_type"].to_numpy() == cell_type
        ct_idx = np.flatnonzero(ct_mask)
        counts = obs["n_counts"].to_numpy()[ct_idx]

        # Condition-independent bin edges from the pooled within-cell-type depths.
        edges = np.quantile(counts, np.linspace(0.0, 1.0, n_bins + 1))
        edges[0], edges[-1] = -np.inf, np.inf
        bin_of = np.digitize(counts, edges[1:-1], right=False)

        conditions = obs["condition"].to_numpy()[ct_idx]
        for b in range(n_bins):
            in_bin = bin_of == b
            per_condition = {
                cond: ct_idx[in_bin & (conditions == cond)]
                for cond in ("control", "TBI")
            }
            n_take = min(len(v) for v in per_condition.values())
            if n_take == 0:
                continue
            for idx in per_condition.values():
                keep.append(rng.choice(idx, size=n_take, replace=False))

    return np.sort(np.concatenate(keep))


def depth_report(obs: pd.DataFrame) -> pd.DataFrame:
    """Median n_counts per cell type and condition, plus the control-TBI gap."""
    med = (
        obs.groupby(["cell_type", "condition"], observed=True)["n_counts"]
        .median()
        .unstack("condition")
    )
    med["abs_gap"] = (med["control"] - med["TBI"]).abs()
    return med


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--in-h5ad",
        required=True,
        help="paper-matched primary AnnData (output of apply_ngene1000_filter.py)",
    )
    ap.add_argument("--out-h5ad", required=True, help="depth-matched AnnData to write")
    ap.add_argument("--n-bins", type=int, default=N_BINS)
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()

    adata = sc.read_h5ad(args.in_h5ad)

    missing = {"cell_type", "condition", "n_counts"} - set(adata.obs.columns)
    if missing:
        raise KeyError(f"input AnnData is missing required obs columns: {sorted(missing)}")

    # The MT corpus is the six major cortical cell types only.
    adata = adata[adata.obs["cell_type"].isin(CELL_TYPES)].copy()
    print(f"before depth matching: {adata.n_obs:,} nuclei")
    print(depth_report(adata.obs).to_string(), end="\n\n")

    keep = depth_match(adata.obs, n_bins=args.n_bins, seed=args.seed)
    matched = adata[keep].copy()

    print(f"after depth matching:  {matched.n_obs:,} nuclei")
    print(matched.obs["condition"].value_counts().to_string())
    print(depth_report(matched.obs).to_string())

    out = Path(args.out_h5ad)
    out.parent.mkdir(parents=True, exist_ok=True)
    matched.write_h5ad(out)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
