#!/usr/bin/env python
"""Eval-only step for the fine-tuned TBI classifier (companion to
12_finetune_tbi_classifier_excitatory.py). Loads the previously-trained
model and runs evaluate_saved_model on the held-out outer test set.

Why separate: the 12_* script trained successfully in job 57711067 (306/306
steps) but its final eval call failed because predict_metadata referenced
'condition' — that column is renamed to numeric 'label' by prepare_data.
Re-running training (~16 min) just to fix the eval call would be wasteful.
"""
from __future__ import annotations
import os, pickle
from pathlib import Path

os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
os.environ.setdefault("WANDB_DISABLED", "true")

from geneformer import Classifier

# Section path configuration (GF_WORK, GF_MODEL_DIR, ...): see ../../gf_paths.py.
import os, sys  # noqa: E401
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import gf_paths  # noqa: E402

BASE = Path(gf_paths.GF_WORK)
# The published run's checkpoint (datestamped 260515 by Geneformer). A rerun of
# 12_finetune_* gets a new datestamp: point GF_FT_MODEL at its ksplit1/checkpoint-<N>.
TRAINED_MODEL = Path(gf_paths.GF_FT_MODEL or gf_paths.FT_PUBLISHED_CHECKPOINT)
ID_CLASS_DICT = BASE / "results/finetune_tbi_excitatory/excit_tbi_clf_id_class_dict.pkl"
TEST_DATASET  = BASE / "results/finetune_tbi_excitatory/excit_tbi_clf_labeled_test.dataset"
HELDOUT_DIR   = BASE / "results/finetune_tbi_excitatory/heldout_eval"
HELDOUT_DIR.mkdir(parents=True, exist_ok=True)

cc = Classifier(
    classifier="cell",
    cell_state_dict={"state_key": "condition", "states": ["TBI", "control"]},
    forward_batch_size=100,
    model_version="V2",
    nproc=8,
    ngpu=1,
)

cc.evaluate_saved_model(
    model_directory=str(TRAINED_MODEL),
    id_class_dict_file=str(ID_CLASS_DICT),
    test_data_file=str(TEST_DATASET),
    output_directory=str(HELDOUT_DIR),
    output_prefix="heldout",
    predict=True,
    predict_metadata=["label", "donor", "age", "sex", "region"],
)

# Load + print the saved metrics
metrics_pkl = HELDOUT_DIR / "heldout_test_metrics_dict.pkl"
with open(metrics_pkl, "rb") as f:
    m = pickle.load(f)

print("=== HELD-OUT TEST METRICS ===")
print(f"Accuracy:  {m.get('acc')}")
print(f"Macro F1:  {m.get('macro_f1')}")
if m.get("all_roc_metrics"):
    print(f"ROC AUC:   {m['all_roc_metrics'].get('all_roc_auc')}")
print(f"Confusion matrix:")
print(m.get("conf_matrix"))

# Pointer for downstream scripts
(BASE / "results/finetune_tbi_excitatory" / "trained_model_path.txt").write_text(str(TRAINED_MODEL) + "\n")
print(f"\nPointer file written: results/finetune_tbi_excitatory/trained_model_path.txt → {TRAINED_MODEL}")
