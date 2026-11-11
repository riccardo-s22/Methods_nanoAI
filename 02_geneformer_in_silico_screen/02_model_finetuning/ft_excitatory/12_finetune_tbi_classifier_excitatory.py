#!/usr/bin/env python
"""Fine-tune Geneformer-V2-104M as a TBI-vs-control cell classifier on
excitatory neurons, using a donor-level held-out split for honest
generalization measurement.

Held-out test donors (never seen during training):
  - control 502 (973 cells, smallest control donor)
  - TBI 6   (275 cells, acute 4 d, model-extrapolation noisy donor)
  - TBI 19  (198 cells, chronic 57 d)

Training pool donors (Geneformer's `validate` does its own random 80/20
inner train/eval split for early stopping; class-stratified):
  - control: 501 (1,478), 529 (1,527)
  - TBI:     1, 2, 3, 7, 8, 10, 11, 16, 20, 21 (1,574 cells total)

Outputs:
  results/finetune_tbi_excitatory/                — prep_data labeled datasets + id_class_dict
  models/Geneformer/fine_tuned_models/<datestamp>_geneformer_cellClassifier_excit_tbi_clf/
                                                  — Geneformer-saved trained model + eval
  results/finetune_tbi_excitatory/heldout_eval/   — confusion matrix, AUC, predictions on held-out donors
"""
from __future__ import annotations

import json
import os
import pickle
import sys
from pathlib import Path

# Reduce noise
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
os.environ.setdefault("WANDB_DISABLED", "true")

from geneformer import Classifier

# Section path configuration (GF_WORK, GF_MODEL_DIR, ...): see ../../gf_paths.py.
import os, sys  # noqa: E401
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import gf_paths  # noqa: E402

BASE = Path(gf_paths.GF_WORK)
TOKENIZED = BASE / "tokenized_paperQC_nGene1000/excitatory_neuron/GSE209552_excitatory_neurons_nGene1000.dataset"
MODEL_ROOT = Path(gf_paths.GF_PRETRAINED)

PREP_OUT_DIR  = BASE / "results/finetune_tbi_excitatory"
FT_PARENT_DIR = Path(gf_paths.GF_FT_MODELS)
HELDOUT_EVAL_DIR = PREP_OUT_DIR / "heldout_eval"

OUTPUT_PREFIX_PREP   = "excit_tbi_clf"
OUTPUT_PREFIX_TRAIN  = "excit_tbi_clf"
OUTPUT_PREFIX_EVAL   = "heldout"

OUTER_TEST_DONORS  = ["502", "6", "19"]
OUTER_TRAIN_DONORS = ["501", "529", "1", "2", "3", "7", "8", "10", "11", "16", "20", "21"]


def main():
    for d in (PREP_OUT_DIR, FT_PARENT_DIR, HELDOUT_EVAL_DIR):
        d.mkdir(parents=True, exist_ok=True)

    print(f"Outer donor split:")
    print(f"  HELD-OUT TEST donors: {OUTER_TEST_DONORS}")
    print(f"  TRAINING POOL donors: {OUTER_TRAIN_DONORS}")

    cc = Classifier(
        classifier="cell",
        cell_state_dict={"state_key": "condition", "states": ["TBI", "control"]},
        training_args={
            "num_train_epochs": 1.0,
            "learning_rate": 5e-5,
            "lr_scheduler_type": "linear",
            "warmup_steps": 500,
            "weight_decay": 0.001,
            "per_device_train_batch_size": 12,
            "per_device_eval_batch_size": 12,
            "logging_steps": 100,
            "save_strategy": "epoch",
            "evaluation_strategy": "epoch",
            "load_best_model_at_end": True,
            "metric_for_best_model": "eval_macro_f1",
            "greater_is_better": True,
            "report_to": "none",
        },
        freeze_layers=2,
        num_crossval_splits=1,
        split_sizes={"train": 0.8, "valid": 0.2, "test": 0.0},
        stratify_splits_col=None,  # HF requires ClassLabel; plain int `label` triggers error.
                                   # Random split is fine: 4,579 train cells / ~2:1 class ratio
                                   # means both splits will retain both classes with high probability.
        forward_batch_size=100,
        model_version="V2",
        nproc=8,
        ngpu=1,
    )

    print("\n=== Step 1/3: prepare_data (donor-level outer split) ===")
    split_id = {
        "attr_key": "donor",
        "train": OUTER_TRAIN_DONORS,
        "test":  OUTER_TEST_DONORS,
    }
    cc.prepare_data(
        input_data_file=str(TOKENIZED),
        output_directory=str(PREP_OUT_DIR),
        output_prefix=OUTPUT_PREFIX_PREP,
        split_id_dict=split_id,
    )

    train_dataset = PREP_OUT_DIR / f"{OUTPUT_PREFIX_PREP}_labeled_train.dataset"
    test_dataset  = PREP_OUT_DIR / f"{OUTPUT_PREFIX_PREP}_labeled_test.dataset"
    id_class_dict = PREP_OUT_DIR / f"{OUTPUT_PREFIX_PREP}_id_class_dict.pkl"

    # Geneformer creates these as directories (HF datasets); confirm existence.
    for p in (train_dataset, test_dataset, id_class_dict):
        if not p.exists():
            sys.exit(f"prepare_data did not produce expected output: {p}")

    print("\n=== Step 2/3: validate (fine-tune classifier on train pool) ===")
    all_metrics_train = cc.validate(
        model_directory=str(MODEL_ROOT),
        prepared_input_data_file=str(train_dataset),
        id_class_dict_file=str(id_class_dict),
        output_directory=str(FT_PARENT_DIR) + "/",
        output_prefix=OUTPUT_PREFIX_TRAIN,
        save_eval_output=True,
        predict_eval=True,
    )
    print("\nInner-validation metrics (random 80/20 split within outer-train pool):")
    print(json.dumps({k: str(v) for k, v in all_metrics_train.items() if k != "conf_matrix"}, indent=2))

    # Locate the actual trained-model directory (datestamped). Pick the most recent.
    matches = sorted(
        [p for p in FT_PARENT_DIR.iterdir()
         if p.is_dir() and OUTPUT_PREFIX_TRAIN in p.name],
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    if not matches:
        sys.exit(f"Could not locate trained-model directory in {FT_PARENT_DIR}")
    trained_dir = matches[0]
    print(f"\nTrained model directory: {trained_dir}")

    # Inside trained_dir there's a ksplit1/ subdir with the checkpoint
    ksplit1 = trained_dir / "ksplit1"
    if ksplit1.exists():
        model_for_eval = ksplit1
    else:
        model_for_eval = trained_dir
    print(f"Loading from: {model_for_eval}")

    print("\n=== Step 3/3: evaluate_saved_model on held-out donors ===")
    cc.evaluate_saved_model(
        model_directory=str(model_for_eval),
        id_class_dict_file=str(id_class_dict),
        test_data_file=str(test_dataset),
        output_directory=str(HELDOUT_EVAL_DIR),
        output_prefix=OUTPUT_PREFIX_EVAL,
        predict=True,
        predict_metadata=["label", "donor", "age", "sex", "region"],  # prepare_data dropped `condition`
    )

    # Load the saved test metrics and summarize
    metrics_pkl = HELDOUT_EVAL_DIR / f"{OUTPUT_PREFIX_EVAL}_test_metrics_dict.pkl"
    if metrics_pkl.exists():
        with open(metrics_pkl, "rb") as f:
            heldout = pickle.load(f)
        print(f"\n=== HELD-OUT TEST METRICS ===")
        print(f"  Accuracy:     {heldout.get('acc')}")
        print(f"  Macro F1:     {heldout.get('macro_f1')}")
        if heldout.get("all_roc_metrics") is not None:
            auc = heldout["all_roc_metrics"].get("all_roc_auc")
            print(f"  ROC AUC:      {auc}")
        if "conf_matrix" in heldout:
            print(f"  Confusion matrix:\n{heldout['conf_matrix']}")

    # Pointer file: record the trained-model path for downstream scripts.
    pointer = PREP_OUT_DIR / "trained_model_path.txt"
    pointer.write_text(str(model_for_eval) + "\n")
    print(f"\nWrote pointer: {pointer}")


if __name__ == "__main__":
    main()
