"""Per-cell-type config for the fine-tuned screen pipeline.

Shared by stages 4/5/6:
  - 05_get_perturb_stats_finetuned.py
  - 07_biological_candidate_filter_finetuned.py
  - 09_pathway_enrichment_finetuned.py

Each entry corresponds to a (cell type, perturbation direction) combination.
The fine-tuned classifier model is symmetric (it learned both axes), so the
same checkpoint can serve both directions; only the perturbation worker's
filter + state choice differs.

Adding a new entry:
  1. Add an entry with tokenized path, run_root, celltype_prefix,
     start_state, goal_state, total_start_cells, cells_per_chunk,
     expected_chunks.
  2. The matching upstream scripts (12/13/perturb worker + Slurm
     launchers) must already exist; this config governs only the
     aggregate/filter/enrichment side.

Note on `total_start_cells`: this is the count of START-state cells (the
ones we perturb), NOT the total. The chunk worker filters to
condition==start_state so chunk indexing assumes start-only cells.
"""
from pathlib import Path

# Section path configuration (GF_WORK, GF_MODEL_DIR, ...): see ../../gf_paths.py.
import os, sys  # noqa: E401
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import gf_paths  # noqa: E402

BASE = Path(gf_paths.GF_WORK)

CELLTYPE_CONFIG: dict[str, dict] = {
    # ----- Excitatory neurons, control → TBI direction (INDUCER screen) -----
    # Perturbing control cells, measuring shift toward TBI centroid.
    # Top hits: TBI-driver candidates (loss/gain of function that induces TBI state).
    "excit": {
        "tokenized": BASE / "tokenized_paperQC_nGene1000/excitatory_neuron/GSE209552_excitatory_neurons_nGene1000.dataset",
        "run_root":  BASE / "results/perturb_neurons_finetuned_control_to_tbi",
        "celltype_prefix": "excit_finetuned",
        "start_state": "control",
        "goal_state":  "TBI",
        "total_start_cells": 3978,
        "cells_per_chunk": 56,
        "expected_chunks": 72,
    },
    # ----- Excitatory neurons, TBI → control direction (RESCUE screen) -----
    # Perturbing TBI cells, measuring shift toward control centroid.
    # Top hits: therapeutic-rescue candidates (loss/gain of function that
    # shifts TBI cells back toward control state).
    # Uses the same fine-tuned model checkpoint as 'excit' — the TBI↔control
    # axis is symmetric, so the state embeddings (control/TBI centroids) are
    # shared (read from the control_to_tbi run_root's state_embs cache).
    "excit_rescue": {
        "tokenized": BASE / "tokenized_paperQC_nGene1000/excitatory_neuron/GSE209552_excitatory_neurons_nGene1000.dataset",
        "run_root":  BASE / "results/perturb_neurons_finetuned_tbi_to_control",
        "celltype_prefix": "excit_finetuned_rescue",
        "start_state": "TBI",
        "goal_state":  "control",
        "total_start_cells": 2047,
        "cells_per_chunk": 56,
        "expected_chunks": 37,    # 2047/56 = 36.6 → 37 (chunks 0-36; chunk 36 has 23 cells)
    },
    # ----- Oligodendrocytes, control → TBI direction (INDUCER screen) -----
    "oligo": {
        "tokenized": BASE / "tokenized_paperQC_nGene1000/oligodendrocyte/GSE209552_oligodendrocytes_nGene1000.dataset",
        "run_root":  BASE / "results/perturb_oligo_finetuned_control_to_tbi",
        "celltype_prefix": "oligo_finetuned",
        "start_state": "control",
        "goal_state":  "TBI",
        "total_start_cells": 763,
        "cells_per_chunk": 56,
        "expected_chunks": 14,
    },
}


def get_config(celltype: str) -> dict:
    if celltype not in CELLTYPE_CONFIG:
        raise SystemExit(
            f"Unknown celltype '{celltype}'. "
            f"Choose from: {sorted(CELLTYPE_CONFIG.keys())}"
        )
    return CELLTYPE_CONFIG[celltype]
