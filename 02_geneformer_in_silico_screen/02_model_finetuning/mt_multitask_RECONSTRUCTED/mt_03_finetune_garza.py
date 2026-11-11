"""Multi-task fine-tuning of Geneformer-V2-104M on the Garza cohort (MT model).

Three classification heads read the final-layer CLS token embedding of a shared
transformer trunk initialized with the pretrained Geneformer weights:

    disease_binary   2 classes   TBI vs control            task weight 0.40
    region           3 classes   frontal / temporal /
                                 fronto-parietal           task weight 0.20
    cell_type        6 classes   excit / inhib / oligo /
                                 OPC / astro / microglia    task weight 0.40

The loss is the task-weighted sum of the per-task cross-entropies; there is no
adversarial, gradient-reversal or decorrelation term.

Hyperparameters (Methods, "Multi-task fine-tuning"):
    max learning rate   2e-5          batch size          12
    scheduler           linear w/     layers frozen       6
                        warmup        classifier dropout  0
    optimizer           AdamW         precision           bf16
    warmup steps        500           epochs              3
    weight decay        0.05          total steps         1,338

Checkpoints are saved every 100 steps and each is evaluated on the HELD-OUT
DONORS. The retained checkpoint has the highest mean macro F1 across the three
heads among checkpoints whose disease head has per-class recall >= 0.20 (which
excludes checkpoints that collapse onto one class).

    Published selection: step 600 of 1,338 (6th of 14 checkpoints),
    mean macro F1 0.721 = (0.732 disease + 0.439 region + 0.992 cell type) / 3.

NOTE ON A KNOWN OPTIMISTIC BIAS: because checkpoint selection used the held-out
donors rather than the validation split, the reported MT metrics are selected on
the same nuclei they are reported for. This differs from the FT model, whose
checkpoint was selected on a validation split. The paper states this; it is
repeated here so the script is not read as implementing a clean holdout.

=============================================================================
Rewritten from the Methods after the original MT code was lost. It documents the
procedure but was not used to generate the published results.
=============================================================================
Caveats:

  1. The TRAINING SEED is UNRECORDED. Neither the Methods nor the surviving
     notes state it, so an exact numerical reproduction of the published
     checkpoint cannot be guaranteed from this script alone. SEED below is a
     placeholder carried over from the data-split step, not a recovered value.

  2. Geneformer's MTLClassifier does not expose a save_steps argument or a
     per-checkpoint held-out evaluation loop. The original pipeline wrapped the
     library with its own loop, writing per-step results under
     `heldout_eval_per_step/step<N>/`; that wrapper is not available.
     `select_checkpoint()` below re-implements the documented selection rule
     over whatever per-step evaluation files are present, rather than pretending
     the library performs the sweep.

Original file name: not recorded (the surviving notes reference only the Slurm
job IDs, 59383802 train + 59383803 eval).
=============================================================================
"""

import argparse
import json
from pathlib import Path

import pandas as pd
from geneformer import MTLClassifier

# Section path configuration: see ../../gf_paths.py. The MT_RECON defaults used
# below are a layout chosen for this repository, not the original MT project's.
import os, sys  # noqa: E401
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import gf_paths  # noqa: E402

# UNRECORDED - see caveat 1 above.
SEED = 42

TASK_COLUMNS = ["disease_binary", "region", "cell_type"]

# Order matches TASK_COLUMNS: disease 0.40, region 0.20, cell type 0.40.
TASK_WEIGHTS = [0.40, 0.20, 0.40]

TOTAL_STEPS = 1338
SAVE_STEPS = 100
WARMUP_STEPS = 500

MANUAL_HYPERPARAMETERS = {
    "learning_rate": 2e-5,
    "warmup_ratio": WARMUP_STEPS / TOTAL_STEPS,  # 500 of 1,338 optimizer steps
    "weight_decay": 0.05,
    "dropout_rate": 0.0,
    "lr_scheduler_type": "linear",
    "use_attention_pooling": False,
    "task_weights": TASK_WEIGHTS,
    "max_layers_to_freeze": 6,
}

MIN_PER_CLASS_RECALL = 0.20


def build_classifier(
    train_path: str,
    val_path: str,
    test_path: str,
    pretrained_path: str,
    model_save_path: str,
    results_dir: str,
) -> MTLClassifier:
    return MTLClassifier(
        task_columns=TASK_COLUMNS,
        train_path=train_path,
        val_path=val_path,
        test_path=test_path,
        pretrained_path=pretrained_path,
        model_save_path=model_save_path,
        results_dir=results_dir,
        batch_size=12,
        epochs=3,
        max_layers_to_freeze={"min": 6, "max": 6},
        use_task_weights=True,
        use_manual_hyperparameters=True,
        manual_hyperparameters=MANUAL_HYPERPARAMETERS,
        seed=SEED,
        tensorboard_log_dir=str(Path(results_dir) / "tblogdir"),
    )


def select_checkpoint(heldout_eval_dir: Path) -> pd.DataFrame:
    """Apply the published checkpoint-selection rule to a per-step sweep.

    Expects one `metrics.json` per checkpoint under
    ``heldout_eval_dir/step<N>/metrics.json``, each holding per-task macro F1 and
    the disease head's per-class recalls, e.g.

        {"macro_f1": {"disease_binary": 0.732, "region": 0.439,
                      "cell_type": 0.992},
         "per_class_recall": {"disease_binary": {"control": 0.970,
                                                 "TBI": 0.621}}}

    Returns the sweep as a table, ordered best-first among eligible checkpoints.
    """
    rows = []
    for metrics_file in sorted(heldout_eval_dir.glob("step*/metrics.json")):
        step = int(metrics_file.parent.name.removeprefix("step"))
        metrics = json.loads(metrics_file.read_text())

        macro_f1 = metrics["macro_f1"]
        recalls = metrics["per_class_recall"]["disease_binary"].values()
        rows.append(
            {
                "step": step,
                "mean_macro_f1": sum(macro_f1[t] for t in TASK_COLUMNS) / len(TASK_COLUMNS),
                **{f"macro_f1_{t}": macro_f1[t] for t in TASK_COLUMNS},
                "min_disease_recall": min(recalls),
                "eligible": min(recalls) >= MIN_PER_CLASS_RECALL,
            }
        )

    if not rows:
        raise FileNotFoundError(f"no step*/metrics.json found under {heldout_eval_dir}")

    sweep = pd.DataFrame(rows).sort_values("step").reset_index(drop=True)
    eligible = sweep[sweep["eligible"]]
    if eligible.empty:
        raise RuntimeError(
            f"no checkpoint reached per-class disease recall >= {MIN_PER_CLASS_RECALL}"
        )

    best = eligible.loc[eligible["mean_macro_f1"].idxmax()]
    print(sweep.to_string(index=False))
    print(
        f"\nretained checkpoint: step {int(best['step'])} of {TOTAL_STEPS} "
        f"(mean macro F1 {best['mean_macro_f1']:.3f}; "
        f"{len(eligible)} of {len(sweep)} checkpoints eligible)"
    )
    return sweep.assign(retained=sweep["step"] == best["step"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--train-path", default=gf_paths.MT_RECON["tokenized_train"],
                    help="tokenized train .dataset")
    ap.add_argument("--val-path", default=gf_paths.MT_RECON["tokenized_valid"],
                    help="tokenized validation .dataset")
    ap.add_argument(
        "--test-path",
        default=gf_paths.MT_RECON["tokenized_test"],
        help="tokenized held-out-donor .dataset (also used for checkpoint selection)",
    )
    ap.add_argument("--pretrained-path", default=gf_paths.GF_PRETRAINED, help="Geneformer-V2-104M")
    ap.add_argument("--model-save-path", default=gf_paths.MT_RECON["model_save"])
    ap.add_argument("--results-dir", default=gf_paths.MT_RECON["finetune_results"])
    ap.add_argument(
        "--select-only",
        action="store_true",
        help="skip training; only apply the checkpoint-selection rule to an existing sweep",
    )
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    heldout_eval_dir = results_dir / "heldout_eval_per_step"

    if not args.select_only:
        mc = build_classifier(
            train_path=args.train_path,
            val_path=args.val_path,
            test_path=args.test_path,
            pretrained_path=args.pretrained_path,
            model_save_path=args.model_save_path,
            results_dir=args.results_dir,
        )
        print(f"fine-tuning: 3 heads {TASK_COLUMNS}, task weights {TASK_WEIGHTS}")
        mc.run_manual_tuning()

        # Checkpoints are written every SAVE_STEPS steps; the per-step held-out
        # evaluation that feeds select_checkpoint() is run here (see caveat 2).
        mc.load_and_evaluate_test_model()

    sweep = select_checkpoint(heldout_eval_dir)
    sweep_out = results_dir / "trajectory_summary.tsv"
    sweep.to_csv(sweep_out, sep="\t", index=False)
    print(f"wrote {sweep_out}")


if __name__ == "__main__":
    main()
