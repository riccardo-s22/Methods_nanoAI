"""Build the MT model's TBI and control state embeddings (goal-state centroids).

Methods ("In silico perturbation"):
    "Start and goal state embeddings were defined as the exact mean
    (EmbExtractor.get_state_embs; summary_stat: exact_mean) of the
    second-to-last-layer CLS token embeddings in each disease state. [...] For
    the MT model, state embeddings were computed from the training-split nuclei
    only (control n=3,819; TBI n=2,868), while perturbations were applied to
    held-out nuclei. The same state embeddings served both directions within
    each model, as only the start and goal states are exchanged."

Two points worth keeping in view when reading this script:

  * Centroids come from the TRAIN split while perturbation runs on the TEST
    split. This departs from Theodoris et al., who draw both from the same pool.
    It is the conservative choice - it removes leakage between the centroid and
    the perturbed nuclei - but it is a departure, not a default.

  * The preprint's 5,000-cell random-subsample cap on centroid construction is
    moot here: both states hold fewer than 5,000 training nuclei, so the exact
    mean is taken over all of them.

=============================================================================
Rewritten from the Methods after the original MT code was lost. It documents the
procedure but was not used to generate the published results.
=============================================================================

Original file name: build_state_embs.py
=============================================================================
"""

import argparse
from pathlib import Path

from geneformer import EmbExtractor

STATE_KEY = "disease_binary"
EMB_LAYER = -1  # second-to-last layer, per Methods
FORWARD_BATCH_SIZE = 16


def state_embs_config(start_state: str, goal_state: str) -> dict:
    """cell_states_to_model for one direction.

    Only start_state and goal_state are exchanged between the rescue and
    induction directions; the centroids themselves are shared.
    """
    return {
        "state_key": STATE_KEY,
        "start_state": start_state,
        "goal_state": goal_state,
        "alt_states": [],
    }


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
        help="tokenized TRAIN-split .dataset of the depth-matched MT corpus",
    )
    ap.add_argument("--output-directory", required=True)
    ap.add_argument("--output-prefix", default="garza_mt_state_embs")
    ap.add_argument("--nproc", type=int, default=8)
    args = ap.parse_args()

    out_dir = Path(args.output_directory)
    out_dir.mkdir(parents=True, exist_ok=True)

    embex = EmbExtractor(
        model_type="MTLCellClassifier",
        emb_mode="cls",
        emb_layer=EMB_LAYER,
        summary_stat="exact_mean",
        forward_batch_size=FORWARD_BATCH_SIZE,
        max_ncells=None,  # both states hold < 5,000 training nuclei
        nproc=args.nproc,
    )

    # TBI -> control supplies both centroids; the induction direction reuses them
    # with start and goal exchanged.
    state_embs = embex.get_state_embs(
        cell_states_to_model=state_embs_config(start_state="TBI", goal_state="control"),
        model_directory=args.model_directory,
        input_data_file=args.input_data_file,
        output_directory=str(out_dir),
        output_prefix=args.output_prefix,
    )

    for state, emb in state_embs.items():
        print(f"{state:>8s} centroid: shape {tuple(emb.shape)}")
    print(f"\nwrote state embeddings to {out_dir}/{args.output_prefix}*")


if __name__ == "__main__":
    main()
