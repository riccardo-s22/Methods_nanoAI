#!/usr/bin/env python
"""Fine-tuned TBI → control RESCUE in-silico perturbation — excitatory neurons.

Mirror of `perturb_cells_finetuned_control_to_tbi.py` but reversed direction:
  • loads the SAME fine-tuned excit TBI classifier (the TBI↔control axis is
    symmetric, so the same model serves both directions)
  • filters TBI cells (not control)
  • measures cosine shift from TBI → control state (start=TBI, goal=control)
  • supports both perturb_type modes via PERTURB_TYPE env var

State embeddings are pre-computed by `13_extract_state_embs_finetuned.py`
(the control-to-TBI side); reuse the same file — centroids are identical.

Chunking: 2,047 TBI cells / 56 per chunk = 37 chunks (array 0-36).
Chunks 0-35 cover 56 cells each (2,016 cells); chunk 36 covers the last 31 cells.

Submit via:
  sbatch scripts/run_perturb_array_finetuned_tbi_to_control_delete.sh
  sbatch scripts/run_perturb_array_finetuned_tbi_to_control_overexpress.sh

Interpretation:
  Top mean_shift > 0 hits = candidate THERAPEUTIC RESCUE genes — perturbations
  (delete or overexpress) that shift TBI cells back toward the control state.
"""
import multiprocessing
if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)

import os
import sys
import time
import io as _io
import pickle
from pathlib import Path

# Section path configuration (GF_WORK, GF_FT_MODEL, ...): see ../../gf_paths.py.
import os, sys  # noqa: E401
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import gf_paths  # noqa: E402

BASE            = Path(gf_paths.GF_WORK)
TOKENIZED       = BASE / "tokenized_paperQC_nGene1000/excitatory_neuron/GSE209552_excitatory_neurons_nGene1000.dataset"
# Reuse the existing fine-tuned model (260515_geneformer_cellClassifier_excit_tbi_clf).
# The pointer file lives in the existing control-to-TBI results dir.
POINTER         = BASE / "results/finetune_tbi_excitatory/trained_model_path.txt"
# Reuse the existing state-embs cache — control and TBI centroids are identical
# regardless of which direction we read them in.
STATE_EMBS_PATH = BASE / "results/perturb_neurons_finetuned_control_to_tbi/state_embs/excitatory_finetuned_state_embs.pkl"
RUN_ROOT        = BASE / "results/perturb_neurons_finetuned_tbi_to_control"

TOTAL_TBI_CELLS = 2047
CELLS_PER_CHUNK = 56     # 37 chunks total (array 0-36); chunk 36 has 31 cells
FORWARD_BATCH   = 30
NPROC           = 1

CELL_STATES = {
    "state_key":   "condition",
    "start_state": "TBI",
    "goal_state":  "control",
    "alt_states":  [],
}


def main():
    perturb_type = os.environ.get("PERTURB_TYPE", "").strip()
    if perturb_type not in {"delete", "overexpress"}:
        sys.exit(
            f"ERROR: PERTURB_TYPE env var must be 'delete' or 'overexpress' (got: {perturb_type!r}).\n"
            "Set it in the Slurm launcher (run_perturb_array_finetuned_tbi_to_control_<mode>.sh)."
        )

    output_dir = RUN_ROOT / f"chunks_{perturb_type}"
    output_dir.mkdir(parents=True, exist_ok=True)

    chunk_id   = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
    start_cell = chunk_id * CELLS_PER_CHUNK
    end_cell   = min(start_cell + CELLS_PER_CHUNK, TOTAL_TBI_CELLS)
    n_cells    = end_cell - start_cell

    if n_cells <= 0:
        print(f"Chunk {chunk_id}: no cells to process (start={start_cell} >= {TOTAL_TBI_CELLS}). Exiting.")
        sys.exit(0)

    output_prefix = f"excit_finetuned_rescue_chunk{chunk_id:03d}_cells{start_cell:04d}-{end_cell - 1:04d}"

    import torch
    print(f"torch {torch.__version__}  |  CUDA: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        p = torch.cuda.get_device_properties(0)
        print(f"GPU: {p.name}  VRAM: {p.total_memory / 1e9:.1f} GB")

    print(f"\nPerturb type: {perturb_type}")
    print(f"Direction:    TBI → control (RESCUE screen)")
    print(f"Chunk:        {chunk_id}  (cells {start_cell}–{end_cell - 1}, n={n_cells})")
    print(f"Output dir:   {output_dir}")
    print(f"Prefix:       {output_prefix}")

    if not POINTER.exists() and not gf_paths.GF_FT_MODEL:
        sys.exit(f"ERROR: fine-tuned model pointer not found: {POINTER}")
    # GF_FT_MODEL if set, else the pointer (re-based onto GF_MODEL_DIR when the
    # stored path is that of another machine); see gf_paths.resolve_model_pointer.
    model_dir = Path(gf_paths.resolve_model_pointer(POINTER))
    if not model_dir.exists():
        sys.exit(f"ERROR: fine-tuned model directory not found: {model_dir}")
    if not STATE_EMBS_PATH.exists():
        sys.exit(f"ERROR: state embeddings not found at {STATE_EMBS_PATH}\nRun 13_extract_state_embs_finetuned.py first.")
    if not TOKENIZED.exists():
        sys.exit(f"ERROR: tokenized dataset not found: {TOKENIZED}")

    # Load state embs to CPU. The file is a pickle.dump of a dict of CUDA
    # tensors; default deserialization tries to restore them to cuda:0 at
    # startup and dies cluster-wide if the assigned GPU has a transient NVML
    # hiccup. Patch torch's per-storage loader to force map_location='cpu'.
    # ISP moves tensors to cuda itself (perturber_utils.py:802).
    import torch.storage as _ts_storage
    _orig_load_from_bytes = _ts_storage._load_from_bytes
    _ts_storage._load_from_bytes = lambda b: torch.load(_io.BytesIO(b), map_location="cpu", weights_only=False)
    try:
        with open(STATE_EMBS_PATH, "rb") as f:
            state_embs_dict = pickle.load(f)
    finally:
        _ts_storage._load_from_bytes = _orig_load_from_bytes
    print(f"State embs:   {list(state_embs_dict.keys())} (loaded from cache, CPU)")
    print(f"Model dir:    {model_dir}")

    from geneformer import InSilicoPerturber

    isp = InSilicoPerturber(
        perturb_type=perturb_type,
        perturb_rank_shift=None,
        genes_to_perturb="all",
        combos=0,
        anchor_gene=None,
        model_type="CellClassifier",      # fine-tuned binary classifier
        num_classes=2,
        emb_mode="cls",
        cell_emb_style="mean_pool",
        filter_data={"condition": ["TBI"]},   # start from TBI cells (rescue direction)
        cell_states_to_model=CELL_STATES,
        state_embs_dict=state_embs_dict,
        max_ncells=None,
        cell_inds_to_perturb={"start": start_cell, "end": end_cell},
        emb_layer=-1,
        forward_batch_size=FORWARD_BATCH,
        nproc=NPROC,
        model_version="V2",
    )

    t0 = time.time()
    isp.perturb_data(
        model_directory=str(model_dir),
        input_data_file=str(TOKENIZED),
        output_directory=str(output_dir),
        output_prefix=output_prefix,
    )
    elapsed = time.time() - t0

    print(f"\nChunk {chunk_id} ({perturb_type}, RESCUE) complete.")
    print(f"  Cells:     {n_cells}")
    print(f"  Runtime:   {elapsed / 60:.1f} min  ({elapsed:.0f} s)")
    print(f"  Sec/cell:  {elapsed / n_cells:.1f}")
    print(f"  Output:    {output_dir}/{output_prefix}_*")


if __name__ == "__main__":
    main()
