"""Donor-level train/validation/test split for the MT model.

Methods ("Multi-task fine-tuning"):
    "Train/test splits were defined at the donor level. One control and one TBI
    donor were held out as the test set, chosen to minimize the largest relative
    difference in cell type proportions between the held-out set and the full
    cohort; within the training donors, a random 80%/20% split (seed 42)
    separated training and validation nuclei."

One donor per class is forced by the cohort: Garza has only 3 control donors, so
holding out two would leave fewer than two controls in training and gut the
disease head. The consequence is carried through the whole MT analysis - the
held-out TBI nuclei all come from a single donor, so Geneformer's donor-
consistency safeguard is undefined (NaN) for every gene in the MT screen and is
neither filtered on nor reported.

=============================================================================
RECONSTRUCTED SCRIPT - NOT THE ORIGINAL CODE THAT PRODUCED THE PUBLISHED RESULTS
=============================================================================
See mt_01_depth_match_garza.py for why the original is unavailable and what this
reconstruction is based on.

Original file name: 01c_donor_split_garza.py
=============================================================================
"""

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

SEED = 42
VALID_FRACTION = 0.20
N_HOLDOUT_PER_CLASS = 1


def celltype_proportions(obs: pd.DataFrame) -> pd.Series:
    return obs["cell_type"].value_counts(normalize=True)


def max_relative_deviation(holdout: pd.DataFrame, full: pd.DataFrame) -> float:
    """Largest relative cell-type proportion difference between holdout and cohort.

    This is the quantity the Methods say the held-out pair was chosen to minimize:
    max over cell types of |p_holdout - p_full| / p_full. A cell type absent from
    the held-out set contributes a deviation of 1.0.
    """
    p_full = celltype_proportions(full)
    p_hold = celltype_proportions(holdout).reindex(p_full.index, fill_value=0.0)
    return float(((p_hold - p_full).abs() / p_full).max())


def choose_holdout_donors(obs: pd.DataFrame) -> tuple[str, str]:
    """Pick the (control, TBI) donor pair minimizing max relative deviation."""
    controls = sorted(obs.loc[obs["condition"] == "control", "donor"].unique())
    tbis = sorted(obs.loc[obs["condition"] == "TBI", "donor"].unique())

    scored = []
    for ctrl, tbi in itertools.product(controls, tbis):
        holdout = obs[obs["donor"].isin([ctrl, tbi])]
        if holdout.empty:
            continue
        scored.append((max_relative_deviation(holdout, obs), ctrl, tbi))

    scored.sort()
    best = scored[0]
    print("candidate held-out donor pairs (lowest max relative deviation first):")
    for dev, ctrl, tbi in scored[:5]:
        print(f"  control {ctrl:>6s} + TBI {tbi:>6s}   max rel. deviation {dev:.3f}")
    return best[1], best[2]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-h5ad", required=True, help="depth-matched MT corpus")
    ap.add_argument("--out-h5ad", required=True)
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()

    adata = sc.read_h5ad(args.in_h5ad)
    obs = adata.obs

    missing = {"cell_type", "condition", "donor"} - set(obs.columns)
    if missing:
        raise KeyError(f"input AnnData is missing required obs columns: {sorted(missing)}")

    ctrl_donor, tbi_donor = choose_holdout_donors(obs)
    print(f"\nheld out: control {ctrl_donor}, TBI {tbi_donor}")

    split = pd.Series("train", index=obs.index, dtype=object)
    test_mask = obs["donor"].isin([ctrl_donor, tbi_donor]).to_numpy()
    split[test_mask] = "test"

    # 80/20 nucleus-level split among the TRAINING donors only.
    rng = np.random.default_rng(args.seed)
    train_idx = np.flatnonzero(~test_mask)
    n_valid = int(round(VALID_FRACTION * len(train_idx)))
    valid_idx = rng.choice(train_idx, size=n_valid, replace=False)
    split.iloc[valid_idx] = "valid"

    adata.obs["split"] = pd.Categorical(split, categories=["train", "valid", "test"])

    print("\nnuclei per split:")
    print(adata.obs["split"].value_counts().reindex(["train", "valid", "test"]).to_string())
    print("\nsplit x condition:")
    print(pd.crosstab(adata.obs["split"], adata.obs["condition"]).to_string())

    out = Path(args.out_h5ad)
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
