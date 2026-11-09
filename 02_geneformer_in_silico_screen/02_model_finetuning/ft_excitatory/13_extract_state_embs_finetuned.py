#!/usr/bin/env python
"""Extract TBI / control state embeddings using the fine-tuned TBI classifier.

Mirrors `extract_state_embs_paperQC_nGene1000_*.py` but loads the fine-tuned
checkpoint (model_type="CellClassifier") instead of the pretrained model.
The resulting embeddings live in the TBI-discriminating space learned during
fine-tuning, so cosine shifts measured against them correspond to movement
along the learned disease axis.

Output:
  results/perturb_neurons_finetuned_control_to_tbi/state_embs/excitatory_finetuned_state_embs.pkl
"""
import multiprocessing
if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)

import pickle
import time
from pathlib import Path

BASE = Path("/work/pi_cfallini_uri_edu/Geneformer_TBI")
TOKENIZED = BASE / "tokenized_paperQC_nGene1000/excitatory_neuron/GSE209552_excitatory_neurons_nGene1000.dataset"
POINTER   = BASE / "results/finetune_tbi_excitatory/trained_model_path.txt"
STATE_DIR = BASE / "results/perturb_neurons_finetuned_control_to_tbi/state_embs"

OUTPUT_PREFIX = "excitatory_finetuned"
FORWARD_BATCH = 16
NPROC = 1

CELL_STATES = {
    "state_key":   "condition",
    "start_state": "control",   # control = start (we'll perturb control cells)
    "goal_state":  "TBI",       # goal = TBI
    "alt_states":  [],
}


def main():
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    cache = STATE_DIR / f"{OUTPUT_PREFIX}_state_embs.pkl"

    if not POINTER.exists():
        raise FileNotFoundError(
            f"Pointer to fine-tuned model not found: {POINTER}. "
            f"Run scripts/12_finetune_tbi_classifier_excitatory.py first."
        )
    MODEL_DIR = Path(POINTER.read_text().strip())
    if not MODEL_DIR.exists():
        raise FileNotFoundError(f"Fine-tuned model directory not found: {MODEL_DIR}")

    import torch
    print(f"torch {torch.__version__}  |  CUDA: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        p = torch.cuda.get_device_properties(0)
        print(f"GPU: {p.name}  VRAM: {p.total_memory / 1e9:.1f} GB")

    print(f"\nTokenized:        {TOKENIZED}")
    print(f"Fine-tuned model: {MODEL_DIR}")
    print(f"Output cache:     {cache}")
    print(f"Cell states:      {CELL_STATES}")

    if cache.exists():
        print(f"\n[skip] State embeddings already cached at {cache}")
        with open(cache, "rb") as f:
            d = pickle.load(f)
        for k, v in d.items():
            print(f"  {k}: shape {v.shape}")
        return

    from geneformer import EmbExtractor

    print("\nRunning EmbExtractor.get_state_embs on fine-tuned model ...")
    t0 = time.time()
    embex = EmbExtractor(
        model_type="CellClassifier",   # fine-tuned binary classifier
        num_classes=2,                 # TBI vs. control
        emb_mode="cls",
        filter_data=None,
        max_ncells=None,
        emb_layer=-1,
        emb_label=None,
        forward_batch_size=FORWARD_BATCH,
        nproc=NPROC,
        summary_stat="exact_mean",
        model_version="V2",
    )
    state_embs_dict = embex.get_state_embs(
        cell_states_to_model=CELL_STATES,
        model_directory=str(MODEL_DIR),
        input_data_file=str(TOKENIZED),
        output_directory=str(STATE_DIR),
        output_prefix=OUTPUT_PREFIX,
    )
    print(f"  Done in {(time.time() - t0)/60:.1f} min")
    for k, v in state_embs_dict.items():
        print(f"    {k}: shape {v.shape}")

    with open(cache, "wb") as f:
        pickle.dump(state_embs_dict, f)
    print(f"\nCached: {cache}")


if __name__ == "__main__":
    main()
