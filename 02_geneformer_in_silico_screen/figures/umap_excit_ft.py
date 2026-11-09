#!/usr/bin/env python
"""UMAP of excitatory-neuron CLS embeddings from the excit-specific fine-tuned
TBI classifier (the "FT" model). Per-cell CLS embeddings are extracted in the
learned TBI-vs-control space, then projected to 2D and colored by condition,
donor, and region.

This is the cell-level view of the disease axis the FT model learned (recall
the FT model is trained ONLY on TBI-vs-control within excit, so condition should
be the dominant organizing direction here, unlike the MT model where cell_type
dominates).

Output: results/umap/excit_ft_umap.png  (+ per-cell embeddings parquet)
"""
import multiprocessing
if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("/work/pi_cfallini_uri_edu/Geneformer_TBI")
TOKENIZED = BASE / "tokenized_paperQC_nGene1000/excitatory_neuron/GSE209552_excitatory_neurons_nGene1000.dataset"
POINTER = BASE / "results/finetune_tbi_excitatory/trained_model_path.txt"
OUT_DIR = BASE / "results/umap"
PREFIX = "excit_ft_cellembs"
EMB_CACHE = OUT_DIR / f"{PREFIX}.parquet"
PNG = OUT_DIR / "excit_ft_umap.png"
LABELS = ["condition", "donor", "region", "sample"]
FORWARD_BATCH = 32


def get_embs():
    if EMB_CACHE.exists():
        print(f"[skip extract] loading cached {EMB_CACHE}")
        return pd.read_parquet(EMB_CACHE)
    MODEL_DIR = Path(POINTER.read_text().strip())
    assert MODEL_DIR.exists(), f"missing model {MODEL_DIR}"
    from geneformer import EmbExtractor
    print(f"extracting CLS embeddings from {MODEL_DIR}")
    embex = EmbExtractor(
        model_type="CellClassifier",
        num_classes=2,
        emb_mode="cls",
        emb_layer=-1,
        emb_label=LABELS,
        max_ncells=None,
        forward_batch_size=FORWARD_BATCH,
        nproc=1,
        model_version="V2",
    )
    df = embex.extract_embs(
        model_directory=str(MODEL_DIR),
        input_data_file=str(TOKENIZED),
        output_directory=str(OUT_DIR),
        output_prefix=PREFIX,
    )
    df.to_parquet(EMB_CACHE)
    print(f"cached embeddings -> {EMB_CACHE}  shape={df.shape}")
    return df


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = get_embs()
    emb_cols = [c for c in df.columns if isinstance(c, (int, np.integer))
                or (isinstance(c, str) and c.isdigit())]
    X = df[emb_cols].to_numpy(dtype=np.float32)
    print(f"embedding matrix: {X.shape}")

    import umap
    reducer = umap.UMAP(n_neighbors=30, min_dist=0.3, metric="cosine",
                        random_state=0)
    XY = reducer.fit_transform(X)
    df["UMAP1"], df["UMAP2"] = XY[:, 0], XY[:, 1]
    df[["UMAP1", "UMAP2"] + [c for c in LABELS if c in df.columns]].to_parquet(
        OUT_DIR / "excit_ft_umap_coords.parquet")

    panels = [c for c in ["condition", "donor", "region"] if c in df.columns]
    fig, axes = plt.subplots(1, len(panels), figsize=(6 * len(panels), 5.5))
    if len(panels) == 1:
        axes = [axes]
    for ax, col in zip(axes, panels):
        cats = sorted(df[col].astype(str).unique())
        cmap = (plt.cm.coolwarm if col == "condition"
                else plt.cm.tab20)
        for i, c in enumerate(cats):
            m = df[col].astype(str) == c
            color = (cmap(0.0 if c == "control" else 1.0) if col == "condition"
                     else cmap(i % 20))
            ax.scatter(df.loc[m, "UMAP1"], df.loc[m, "UMAP2"], s=4, alpha=0.5,
                       linewidths=0, color=color, label=f"{c} (n={int(m.sum())})")
        ax.set_title(f"FT excit CLS embedding — by {col}")
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
        ax.legend(markerscale=3, fontsize=7, loc="best", framealpha=0.7)
    fig.suptitle("Excitatory-neuron fine-tuned (TBI-vs-control) model — "
                 f"per-cell CLS embedding UMAP (n={len(df):,})", fontsize=12)
    fig.tight_layout()
    fig.savefig(PNG, dpi=150, bbox_inches="tight")
    print(f"wrote {PNG}")


if __name__ == "__main__":
    main()
