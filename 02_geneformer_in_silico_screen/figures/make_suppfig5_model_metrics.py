#!/usr/bin/env python
"""
Supplementary Figure 5 -- classification performance of the two Geneformer models.

Panels written (editable text, Arial):
  suppfig5_metrics__table.svg   summary table: accuracy, macro F1, ROC AUC
  suppfig5_metrics__cm_FT.svg   FT confusion matrix
  suppfig5_metrics__cm_MT.svg   MT confusion matrix
  suppfig5_metrics__full.svg    the three side by side

DATA PROVENANCE -- every number below is transcribed from a logged evaluation,
nothing is recomputed from a model and nothing is traced off a figure.

  FT (fine-tuned excitatory-neuron TBI-vs-control classifier)
      <project root>/logs/finetune_eval_57711262.out   (job 57711262,
      uri-gpu012, 15 May 2026) section "=== HELD-OUT TEST METRICS ===".
      Design is a donor-level outer split, from finetune_tbi_excit_57711067.out:
          HELD-OUT TEST donors : 502, 6, 19
          TRAINING POOL donors : 501, 529, 1, 2, 3, 7, 8, 10, 11, 16, 20, 21
      so the reported metrics are on nuclei from donors the model never saw.
      The inner-validation numbers (random 80/20 within the training pool,
      acc 0.9607 / macro F1 0.9561 / AUC 0.9878) are NOT plotted here: they are
      the model-selection split, not a performance estimate.

  MT (Garza multi-task model, "MT v1 depth-matched", checkpoint-600)
      The primary eval log is GONE: the project tree
      <cluster work area>/Geneformer_Garza/ no longer exists on Unity, the
      /snapshots/work backups reach back only to 2026-11-05 (already after the
      deletion), and Slurm accounting confirms only that the eval job ran --
          sacct -j 59383803 -> garza_mt_eval  COMPLETED  2026-05-28T13:47:51
                               WorkDir <cluster work area>/Geneformer_Garza/scripts
      The model checkpoint is deleted with it, so the evaluation cannot be re-run.

      The values below are the ones that survive, in two independent Claude-Code
      project-memory records written at the time, on Unity (the coding
      assistant's project-memory folder for the Geneformer_TBI project):
          project_geneformer_garza_mt_v1.md
            "disease_binary 0.732 - per-class recall control 0.970, TBI 0.621"
            "cell_type 0.992", "region 0.439", mean macro-F1 0.721
            eval corpus 9,428 depth-matched nuclei (4,714 control / 4,714 TBI);
            holdout 29% of corpus, donor-level: 1 control donor (502) + 3 TBI donors
          project_isp_shift_geometry_mt_vs_ft.md
            "cell_type recall 0.99 vs TBI recall 0.62"

      ACCURACY, ROC AUC AND THE CONFUSION-MATRIX COUNTS WERE NEVER RECORDED and are
      not recoverable; they are left blank rather than reconstructed. The two recalls
      do fully determine the ROW-NORMALISED matrix, so the MT panel is drawn as
      percentages with the counts marked absent.

      *** DISCREPANCY WITH THE CURRENT FIGURE 2a LEGEND ***
      The legend states condition recall "0.972 for controls and 0.833 for TBI".
      The surviving record of this model gives 0.970 / 0.621. The control figures
      agree; the TBI figure does not, and the two cannot come from one evaluation:
      recalls of 0.972/0.833 imply a macro-F1 of ~0.902 at any class ratio, against
      the recorded 0.732. (Solving the recorded macro-F1 for the holdout class ratio
      gives N_TBI/N_control = 2.07, consistent with the 1-control/3-TBI donor split.)
      project_geneformer_garza_mt_v1.md also warns that TBI recall ~0.97 is what the
      NON-depth-matched arms produced, and that depth-matching drops it to ~0.62 --
      "the 35-point drop is the magnitude of the depth shortcut that had been
      masquerading as disease signal". 0.833 is consistent with neither. Resolve
      which model the legend describes before this panel is published.

WHAT IS DERIVED, AND HOW (marked with a dagger in the table)
      For a 2-class problem the per-class recalls rc, rt plus the macro-F1 leave
      exactly one free parameter -- the holdout class ratio r = N_TBI / N_control --
      so r, and with it the accuracy, are IDENTIFIED:
          macro-F1(rc, rt, r) = 0.732   ->   r = 2.067
          accuracy = (rc + rt*r) / (1 + r) = 0.735
      Running the whole solve over the rounding intervals of all three recorded
      values (each is given to 3 dp, so each is +-0.0005) moves these only to
          r        in [2.046, 2.088]
          accuracy in [0.7343, 0.7353]      -> 0.735 is safe to 3 dp
      r = 2.07 is an independent sanity check on the record: the holdout was
      1 control donor against 3 TBI donors, so more TBI nuclei than control is
      exactly what the donor split predicts.
      Balanced accuracy needs no r at all: (rc + rt)/2 = 0.796, exact.

      ROC AUC IS NOT IDENTIFIED and is not reported as a number. A confusion
      matrix is one operating point (FPR 0.030, TPR 0.621); the AUC depends on the
      whole score distribution. Any monotone ROC through that point has AUC in
      [0.602, 0.989]; assuming a concave (proper) ROC gives only the lower bound
      AUC >= 0.795, which is just the balanced accuracy again. A bound is what the
      table shows.

      The implied counts (N_control ~892, N_TBI ~1843, from 29% of the 9,428-nucleus
      corpus) are shown with "~" in the MT panel and nowhere else: "29%" is itself
      rounded, which alone moves them to [876, 907] and [1811, 1874]. They are an
      order-of-magnitude aid, not data.

TO COMPLETE THE FIGURE: if the MT eval output is recovered, fill accuracy / roc_auc /
cm for MT below and re-run; panels that have data are drawn identically either way.
Recomputing the derivation: scratchpad mt_derive.py.

Output folder in this repository: <results>/02_geneformer_in_silico_screen/figures/
(originally next to this script). The script reads no input file.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Section path configuration (GF_WORK, results folder, shipped data): see ../gf_paths.py.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gf_paths  # noqa: E402

STEM = os.path.join(gf_paths.gf_results_dir("figures"), "suppfig5_metrics")

# ----------------------------------------------------------------- the numbers
# cm      : rows = true class, cols = predicted class, in CLASSES order (counts).
# recall  : fallback when counts were never recorded -- fixes the row-normalised
#           matrix exactly, which is all the panel's colour and % labels need.
MODELS = {
    "FT": dict(
        title="excitatory-neuron classifier",
        sub="held-out donors 502, 6, 19  (n = 1,446 nuclei)",
        color="#2a78d6",
        accuracy=0.9446749654218534,
        macro_f1=0.9340672869326925,
        roc_auc=0.9782390940162398,
        cm=np.array([[973, 0],
                     [80, 393]]),
        recall=None,
        auc_bound=None,
        approx_n=None,
        derived=set(),
    ),
    "MT": dict(
        title="multi-task model, condition head",
        sub="depth-matched, checkpoint-600  (counts not recorded)",
        color="#eb6834",
        accuracy=0.7348,            # DERIVED, see docstring
        macro_f1=0.732,
        roc_auc=None,               # not identified by one operating point
        auc_bound=0.795,            # DERIVED lower bound (concave ROC)
        cm=None,
        recall=[0.970, 0.621],
        approx_n=[892, 1843],       # DERIVED, "~" only
        # balanced accuracy is NOT derived: it is exactly (rc + rt)/2 of the two
        # recorded recalls, so it needs no reconstruction and carries no dagger.
        derived={"accuracy", "auc_bound"},
    ),
}
CLASSES = ["control", "TBI"]

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
plt.rcParams.update({
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.facecolor": "white",
})

NA = "n.d."


def _white_on(frac):
    """Label ink for a cell: white only once the fill is dark enough to need it."""
    return "white" if frac > 0.55 else INK


def _heads(ax, key):
    """Panel title + subtitle, stacked above the axes so neither can collide.
    The model tag carries the model's colour, matching the table's column heads;
    the words beside it name the model, so identity never rests on colour alone."""
    m = MODELS[key]
    t = ax.text(0, 1.20, key, transform=ax.transAxes, fontsize=9,
                color=m["color"], weight="bold", va="baseline")
    ax.annotate(" " + m["title"], xycoords=t, xy=(1, 0), fontsize=9,
                color=INK, va="baseline")
    ax.text(0, 1.07, m["sub"], transform=ax.transAxes, fontsize=7.5,
            color=MUTED, va="baseline")


def confusion_panel(ax, key):
    """2x2 confusion matrix. Fill = row-normalised recall (sequential, one hue,
    light->dark), so colour encodes the rate while the printed count stays absolute."""
    m = MODELS[key]
    _heads(ax, key)

    cm = m["cm"]
    n = len(CLASSES)
    # right margin inside the axes holds the per-row n, so it cannot reach the
    # neighbouring panel
    ax.set_xlim(-0.56, n - 0.5 + 0.70); ax.set_ylim(n - 0.44, -0.56)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(CLASSES, fontsize=8)
    ax.set_yticklabels(CLASSES, fontsize=8)
    ax.set_xlabel("predicted", fontsize=8, color=INK2, labelpad=4)
    ax.set_ylabel("true", fontsize=8, color=INK2, labelpad=4)
    ax.tick_params(length=0, pad=3)
    for s in ax.spines.values():
        s.set_visible(False)

    if cm is None and m["recall"] is None:
        ax.add_patch(plt.Rectangle((-0.46, -0.46), n - 0.08, n - 0.08,
                                   facecolor="#f7f6f3", edgecolor=GRID, lw=0.8,
                                   hatch="////", zorder=1))
        ax.text(0.5 * (n - 1), 0.5 * (n - 1), "not available",
                ha="center", va="center", fontsize=8.5, color=MUTED, zorder=2)
        return

    # row-normalised matrix: from counts when we have them, else from the recalls
    if cm is not None:
        tot = cm.sum(1, keepdims=True)
        rown = cm / tot
    else:
        r = np.asarray(m["recall"], float)
        rown = np.column_stack([r, 1 - r])          # 2 classes: off-diagonal = 1 - recall
        rown[1] = rown[1][::-1]                      # row "TBI": [1-recall, recall]
        tot = None

    base = matplotlib.colors.to_rgb(m["color"])
    for i in range(n):
        for j in range(n):
            frac = rown[i, j]
            # light -> dark on a single hue: blend the hue toward white
            fill = tuple(1 - frac * (1 - c) for c in base)
            # faint ring so a 0-count (white) cell is still a visible cell
            ax.add_patch(plt.Rectangle((j - 0.46, i - 0.46), 0.92, 0.92,
                                       facecolor=fill, edgecolor=GRID, lw=0.8, zorder=1))
            if tot is not None:
                ax.text(j, i - 0.09, "%d" % cm[i, j], ha="center", va="center",
                        fontsize=12, color=_white_on(frac), zorder=2)
                ax.text(j, i + 0.20, "%.1f%%" % (100 * frac), ha="center", va="center",
                        fontsize=8, color=_white_on(frac), zorder=2)
            else:
                ax.text(j, i, "%.1f%%" % (100 * frac), ha="center", va="center",
                        fontsize=12, color=_white_on(frac), zorder=2)
    for i in range(n):
        if tot is not None:
            lab = "n = %d" % tot[i, 0]
        elif m["approx_n"] is not None:
            lab = "n ≈ %d†" % m["approx_n"][i]      # derived, never measured
        else:
            lab = "n.r."
        ax.text(n - 0.44, i, lab, ha="left", va="center", fontsize=7.5, color=MUTED)


# Rows both models can fill WITHOUT a reconstruction: macro F1 and the two recalls
# are recorded for each, and balanced accuracy is an exact function of the two
# recalls, so it is measured-equivalent rather than derived. Accuracy and ROC AUC
# are excluded -- for MT the first is reconstructed and the second only bounded, so
# putting them beside FT's measured values would invite a false comparison.
ROWS_SHARED = [("Balanced accuracy", "bal_acc"), ("Macro F1", "macro_f1"),
               ("Recall, control", "rec_c"), ("Recall, TBI", "rec_t")]
ROWS_ALL = [("Accuracy", "accuracy"), ("Balanced accuracy", "bal_acc"),
            ("Macro F1", "macro_f1"), ("ROC AUC", "roc_auc"),
            ("Recall, control", "rec_c"), ("Recall, TBI", "rec_t")]


def table_panel(ax, rows=None):
    """Summary table. Values are text, so they stay selectable in the SVG."""
    ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    keys = list(MODELS)
    rows = ROWS_SHARED if rows is None else rows
    shared = rows is ROWS_SHARED
    x0, x1, xs = 0.015, 0.985, [0.60, 0.86]

    def rec(m, i):
        if m["cm"] is not None:
            return m["cm"][i, i] / m["cm"][i].sum()
        return None if m["recall"] is None else m["recall"][i]

    def val(k, fld):
        """-> (text, is_derived). None text means 'not available'."""
        m = MODELS[k]
        dag = fld in m["derived"]
        if fld == "rec_c":
            v = rec(m, 0)
        elif fld == "rec_t":
            v = rec(m, 1)
        elif fld == "bal_acc":
            a, b = rec(m, 0), rec(m, 1)
            v = None if a is None or b is None else 0.5 * (a + b)
        elif fld == "roc_auc":
            if m["roc_auc"] is not None:
                return "%.3f" % m["roc_auc"], False
            if m["auc_bound"] is not None:
                return "≥ %.3f†" % m["auc_bound"], True
            return None, False
        else:
            v = m[fld]
        if v is None:
            return None, False
        return ("%.3f†" % v) if dag else ("%.3f" % v), dag

    ax.text(x0, 1.20, "Classification performance", fontsize=9, color=INK,
            va="baseline")
    for k, x in zip(keys, xs):
        ax.text(x, 0.95, k, fontsize=9, color=MODELS[k]["color"],
                ha="center", va="baseline", weight="bold")
    ax.plot([x0, x1], [0.905] * 2, color=INK2, lw=0.9)

    dy = 0.128 if len(rows) > 4 else 0.165
    ytop = 0.815 if len(rows) > 4 else 0.78
    for r, (lab, fld) in enumerate(rows):
        y = ytop - r * dy
        ax.text(x0, y, lab, fontsize=8.5, color=INK2, va="center")
        for k, x in zip(keys, xs):
            txt, dag = val(k, fld)
            ax.text(x, y, NA if txt is None else txt, fontsize=10,
                    color=MUTED if (txt is None or dag) else INK,
                    ha="center", va="center")
        if r < len(rows) - 1:
            ax.plot([x0, x1], [y - dy / 2] * 2, color=GRID, lw=0.7)
    ax.plot([x0, x1], [ytop - (len(rows) - 0.5) * dy] * 2, color=INK2, lw=0.9)

    notes = ((
        "Metrics recorded for both models. Balanced accuracy is the mean of",
        "the two recalls. Accuracy and ROC AUC are omitted: for MT they were",
        "never recorded (see script docstring).",
    ) if shared else (
        "† derived from the recorded recalls and macro F1, not measured.",
        "   AUC is bounded, not identified: one operating point cannot fix it.",
        "   The MT eval log and checkpoint were deleted (see script docstring).",
    ))
    for i, t in enumerate(notes):
        ax.text(x0, -0.005 - 0.075 * i, t, fontsize=6.8, color=MUTED, va="center")


def save(fig, name):
    for ext in ("svg", "png"):
        fig.savefig("%s__%s.%s" % (STEM, name, ext))
    plt.close(fig)
    print("wrote %s__%s.svg / .png" % (STEM, name))


# A common top margin for every panel, so the three sit on one baseline when the
# composite is assembled and when the panels are placed individually in the deck.
RECT = (0.16, 0.17, 0.80, 0.66)

# ------------------------------------------------------------------- individual
for k in MODELS:
    f = plt.figure(figsize=(2.9, 2.7))
    confusion_panel(f.add_axes(RECT), k)
    save(f, "cm_%s" % k)

# "table"      -> shared metrics only (the one meant for the figure)
# "table_full" -> every metric, MT's reconstructions daggered (for the response
#                 to reviewers / internal record, not for the panel)
for nm, rws in (("table", ROWS_SHARED), ("table_full", ROWS_ALL)):
    f = plt.figure(figsize=(3.5, 3.0))
    table_panel(f.add_axes((0.04, 0.17, 0.92, 0.66)), rws)
    save(f, nm)

# ------------------------------------------------------------------------ full
for nm, rws in (("full", ROWS_SHARED), ("full_all_metrics", ROWS_ALL)):
    fig = plt.figure(figsize=(9.2, 2.7))
    w, gap, left = 0.235, 0.075, 0.055
    confusion_panel(fig.add_axes((left, 0.17, w, 0.66)), "FT")
    confusion_panel(fig.add_axes((left + w + gap, 0.17, w, 0.66)), "MT")
    table_panel(fig.add_axes((left + 2 * (w + gap), 0.17, w + 0.03, 0.66)), rws)
    save(fig, nm)
