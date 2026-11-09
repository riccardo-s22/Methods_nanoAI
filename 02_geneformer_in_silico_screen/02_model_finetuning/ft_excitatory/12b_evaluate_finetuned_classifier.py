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

BASE = Path("/work/pi_cfallini_uri_edu/Geneformer_TBI")
TRAINED_MODEL = BASE / "models/Geneformer/fine_tuned_models_tbi/260515_geneformer_cellClassifier_excit_tbi_clf/ksplit1/checkpoint-306"
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
