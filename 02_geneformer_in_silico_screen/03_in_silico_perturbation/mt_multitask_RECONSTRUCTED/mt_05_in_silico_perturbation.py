"""In silico perturbation with the MT model - all four arms.

Methods ("In silico perturbation", Table 1 and the MT paragraph):

    arm                 start cells   goal state   perturbation   positive hit means
    rescue / delete        TBI          control      delete       losing the gene moves TBI cells toward control
    rescue / overexpress   TBI          control      overexpress  activating the gene moves TBI cells toward control
    induction / delete     control      TBI          delete       losing the gene drives control cells toward TBI
    induction / overexpress control      TBI          overexpress  activating the gene drives control cells toward TBI

    "The same four arms were run with the MT model, using the TBI and control
    state embeddings derived from its disease state task across all cell types.
    All start-state nuclei of the held-out donors were perturbed. The
    TBI-to-control rescue arms comprised 1,846 held-out TBI nuclei, all from a
    single donor because the split held out one donor per condition. Shifts were
    stratified by annotated cell type, and per-gene statistics were computed
    separately for each cell type with the same test, safeguards and thresholds
    as for the FT model; excitatory neurons are reported in the main text."

Unlike the FT arms, the MT start-state nuclei were NOT subsampled: every
held-out nucleus of the start state is perturbed.

=============================================================================
Rewritten from the Methods after the original MT code was lost. It documents the
procedure but was not used to generate the published results.
=============================================================================

UNRECORDED: the chunk size used to fit the MT arms inside the cluster wall-time
limit is not stated for the MT model (the Methods give 56 nuclei per chunk for
the FT model only). --chunk-size below defaults to the FT value so the script is
runnable; it is an assumption.

Original file name: not recorded.
=============================================================================
"""

import argparse
import pickle
from pathlib import Path

from geneformer import InSilicoPerturber

STATE_KEY = "disease_binary"
EMB_LAYER = -1
FORWARD_BATCH_SIZE = 30
CHUNK_SIZE = 56  # UNRECORDED for the MT model; FT value, see header

ARMS = {
    "rescue_delete": ("TBI", "control", "delete"),
    "rescue_overexpress": ("TBI", "control", "overexpress"),
    "induction_delete": ("control", "TBI", "delete"),
    "induction_overexpress": ("control", "TBI", "overexpress"),
}


def load_state_embs(path: Path) -> dict:
    with open(path, "rb") as handle:
        return pickle.load(handle)


def run_arm(
    arm: str,
    model_directory: str,
    input_data_file: str,
    output_directory: Path,
    state_embs_dict: dict,
    chunk_index: int | None,
    chunk_size: int,
    nproc: int,
) -> None:
    start_state, goal_state, perturb_type = ARMS[arm]

    # Slurm array tasks each take one fixed chunk of start-state nuclei; shifts
    # are written per chunk and pooled before statistics.
    if chunk_index is None:
        cell_inds_to_perturb = "all"
        out_prefix = f"garza_mt_{arm}"
    else:
        cell_inds_to_perturb = {
            "start": chunk_index * chunk_size,
            "end": (chunk_index + 1) * chunk_size,
        }
        out_prefix = f"garza_mt_{arm}_chunk{chunk_index:03d}"

    isp = InSilicoPerturber(
        perturb_type=perturb_type,
        genes_to_perturb="all",
        combos=0,
        model_type="MTLCellClassifier",
        emb_mode="cls",
        cell_emb_style="mean_pool",
        emb_layer=EMB_LAYER,
        filter_data={STATE_KEY: [start_state]},
        cell_states_to_model={
            "state_key": STATE_KEY,
            "start_state": start_state,
            "goal_state": goal_state,
            "alt_states": [],
        },
        state_embs_dict=state_embs_dict,
        max_ncells=None,  # every held-out start-state nucleus is perturbed
        cell_inds_to_perturb=cell_inds_to_perturb,
        forward_batch_size=FORWARD_BATCH_SIZE,
        nproc=nproc,
    )

    arm_dir = output_directory / arm
    arm_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{arm}] {start_state} -> {goal_state}, {perturb_type}, cells={cell_inds_to_perturb}")
    isp.perturb_data(
        model_directory=model_directory,
        input_data_file=input_data_file,
        output_directory=str(arm_dir),
        output_prefix=out_prefix,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", required=True, choices=sorted(ARMS))
    ap.add_argument(
        "--model-directory",
        required=True,
        help="retained MT checkpoint (published: .../ksplit1/checkpoint-600)",
    )
    ap.add_argument(
        "--input-data-file",
        required=True,
        help="tokenized HELD-OUT-DONOR .dataset (the perturbed nuclei)",
    )
    ap.add_argument(
        "--state-embs",
        required=True,
        help="pickled centroids written by mt_04_build_state_embs.py",
    )
    ap.add_argument("--output-directory", required=True)
    ap.add_argument(
        "--chunk-index",
        type=int,
        default=None,
        help="Slurm array index; omit to perturb all start-state nuclei in one job",
    )
    ap.add_argument("--chunk-size", type=int, default=CHUNK_SIZE)
    ap.add_argument("--nproc", type=int, default=8)
    args = ap.parse_args()

    run_arm(
        arm=args.arm,
        model_directory=args.model_directory,
        input_data_file=args.input_data_file,
        output_directory=Path(args.output_directory),
        state_embs_dict=load_state_embs(Path(args.state_embs)),
        chunk_index=args.chunk_index,
        chunk_size=args.chunk_size,
        nproc=args.nproc,
    )


if __name__ == "__main__":
    main()
