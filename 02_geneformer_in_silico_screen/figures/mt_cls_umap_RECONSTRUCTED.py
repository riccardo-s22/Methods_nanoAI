"""Figure 2a - UMAP of MT-model CLS embeddings.

Extracts the second-to-last-layer CLS token embedding of every nucleus in the
depth-matched MT corpus and reduces it to two dimensions, colored by disease
state and outlined by annotated cell type.

Methods ("Cell embeddings and UMAP"):
    "Cell embeddings were extracted as the second-to-last-layer CLS token
    embedding (768 embedding dimensions) of each nucleus by forward pass through
    the respective fine-tuned model (EmbExtractor; model_type: CellClassifier;
    emb_mode: cls; emb_layer: -1; forward batch size: 32; max_ncells: None).
    Cell embeddings from the MT model were extracted for the 9,428 disease
    state-balanced, depth-matched nuclei of the MT corpus [...] Embeddings were
    reduced to two dimensions with UMAP (umap-learn; n_neighbors: 30;
    min_dist: 0.3; metric: cosine; random_state: 0) and colored by disease state
    and outlined by cell type."

Published panel composition (Figure 2 legend): 9,428 condition-balanced nuclei
(4,714 TBI, 4,714 control) - excitatory neurons 3,372; inhibitory 3,078;
oligodendrocytes 1,418; OPCs 716; astrocytes 492; microglia/macrophages 352.

Note that `emb_layer: -1` is the SECOND-TO-LAST layer in Geneformer's
convention, not the last. The classification heads read the final layer; the
embeddings plotted here come from one layer earlier.

All 9,428 nuclei are embedded - this is the whole depth-matched corpus, not a
sample drawn for the figure.

=============================================================================
RECONSTRUCTED SCRIPT - NOT THE ORIGINAL CODE THAT PRODUCED THE PUBLISHED FIGURE
=============================================================================
The MT pipeline lived in two cluster project directories (Geneformer_Garza and
Geneformer_JointPG) that were deleted before this repository was assembled; no
copy survives on the attached drives, in the cluster home directory, or in the
three retained /work snapshots, all of which postdate the deletion. The MT
checkpoint itself is gone with them, so this script cannot currently be run to
completion even with the data in hand.

Rewritten from the published Methods and figure legend and against the
Geneformer API as vendored under models/Geneformer. The sibling
`umap_excit_ft.py` in this directory is the ORIGINAL script for the equivalent
FT panel (Supplementary Figure 5b) and was used as the template for the plotting
conventions here.

UNRECORDED: the exact styling of the published panel (marker sizes, palette,
"outlined by cell type" rendering). What follows is a faithful implementation of
the described content, not a pixel reproduction.

Original file name: not recorded.
=============================================================================
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

EMB_LAYER = -1  # second-to-last layer
FORWARD_BATCH_SIZE = 32

UMAP_KWARGS = dict(n_neighbors=30, min_dist=0.3, metric="cosine", random_state=0)

CONDITION_COLORS = {"control": "#4878a8", "TBI": "#c84848"}

CELL_TYPE_ORDER = [
    "excitatory_neuron",
    "inhibitory_neuron",
    "oligodendrocyte",
    "OPC",
    "astrocyte",
    "microglia_macrophage",
]


def extract_cls_embeddings(model_directory: str, input_data_file: str, out_dir: Path,
                           nproc: int) -> pd.DataFrame:
    """Per-nucleus CLS embeddings with their condition and cell-type labels."""
    from geneformer import EmbExtractor

    embex = EmbExtractor(
        model_type="MTLCellClassifier",
        emb_mode="cls",
        emb_layer=EMB_LAYER,
        emb_label=["disease_binary", "cell_type", "donor"],
        forward_batch_size=FORWARD_BATCH_SIZE,
        max_ncells=None,  # all 9,428 nuclei
        nproc=nproc,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    return embex.extract_embs(
        model_directory=model_directory,
        input_data_file=input_data_file,
        output_directory=str(out_dir),
        output_prefix="garza_mt_cls_embs",
    )


def plot_umap(coords: np.ndarray, labels: pd.DataFrame, out_path: Path) -> None:
    """Colored by disease state, outlined by cell type."""
    fig, ax = plt.subplots(figsize=(6.2, 5.6))

    for condition, color in CONDITION_COLORS.items():
        mask = (labels["disease_binary"] == condition).to_numpy()
        ax.scatter(
            coords[mask, 0], coords[mask, 1],
            s=3, c=color, alpha=0.55, linewidths=0, label=f"{condition} (n={mask.sum():,})",
        )

    # Cell-type outlines: the convex-hull-free version of "outlined by cell
    # type" - a labelled contour at each cell type's centroid.
    present = [ct for ct in CELL_TYPE_ORDER if (labels["cell_type"] == ct).any()]
    for cell_type in present:
        mask = (labels["cell_type"] == cell_type).to_numpy()
        cx, cy = coords[mask, 0].mean(), coords[mask, 1].mean()
        ax.annotate(
            f"{cell_type.replace('_', ' ')}\nn={mask.sum():,}",
            (cx, cy), ha="center", va="center", fontsize=7,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#555555", alpha=0.75, lw=0.6),
        )

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("MT model, second-to-last-layer CLS embeddings", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8, markerscale=3, loc="best")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.tight_layout()
    for suffix in (".png", ".svg", ".pdf"):
        fig.savefig(out_path.with_suffix(suffix), dpi=400, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--model-directory",
        required=True,
        help="retained MT checkpoint (published: .../ksplit1/checkpoint-600)",
    )
    ap.add_argument(
        "--input-data-file",
        required=True,
        help="tokenized depth-matched MT corpus (all 9,428 nuclei)",
    )
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--nproc", type=int, default=8)
    args = ap.parse_args()

    import umap

    out_dir = Path(args.out_dir)
    embs = extract_cls_embeddings(
        args.model_directory, args.input_data_file, out_dir, args.nproc
    )

    label_cols = ["disease_binary", "cell_type", "donor"]
    labels = embs[label_cols].reset_index(drop=True)
    matrix = embs.drop(columns=label_cols).to_numpy(dtype=np.float32)
    print(f"embeddings: {matrix.shape[0]:,} nuclei x {matrix.shape[1]} dimensions")
    print(labels["disease_binary"].value_counts().to_string())
    print(labels["cell_type"].value_counts().to_string())

    coords = umap.UMAP(**UMAP_KWARGS).fit_transform(matrix)

    pd.DataFrame(
        {"umap1": coords[:, 0], "umap2": coords[:, 1], **{c: labels[c] for c in label_cols}}
    ).to_csv(out_dir / "garza_mt_cls_umap_coords.tsv.gz", sep="\t", index=False)

    plot_umap(coords, labels, out_dir / "fig2a_mt_cls_umap")
    print(f"wrote {out_dir}/fig2a_mt_cls_umap.{{png,svg,pdf}}")


if __name__ == "__main__":
    main()
