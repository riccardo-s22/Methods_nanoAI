#!/usr/bin/env python
"""
Standalone, self-contained SVG of Figure 2b -- cross-model consensus compounds.

WHY THIS EXISTS, AND HOW IT DIFFERS FROM excit_drug_repurpose_figure__panel_b.svg
  That file is the faithful panel B of the 3-panel composite, regenerated from the
  original script unchanged (via ../make_svg.py). In main_figures.pptx the composite's
  panel B is placed as a CROPPED RASTER (srcRect l=13.0% t=48.1% r=49.4% b=4.6%), which
  cuts off the panel title, the drug labels and the x-axis label; those were then re-added
  by hand as 13 separate PowerPoint text boxes.
  This script rebuilds panel B as ONE self-contained vector panel carrying all of its own
  text, so the cropped raster and all 13 overlay text boxes can be deleted and replaced
  by a single SVG.

DATA: the same file the published figure used --
  <project root>/results/drug_repurpose/crossmodel_consensus.tsv
  In this repository it is looked up below
  <results>/02_geneformer_in_silico_screen/drug_repurpose/ (a rerun of
  excit_drug_crossmodel.py), then below $GF_WORK/results/drug_repurpose/, then in
  data/02_geneformer_in_silico_screen/drug_signature_enrichment/ (shipped).
  (the composite PNG in main_figures.pptx is md5 c2d38bae5246680a4a76838853ab2333,
   identical to results/drug_repurpose/excit_drug_repurpose_figure.png)
Nothing is traced off the raster and no value is retyped.

PRESENTATION follows the manuscript caption of Figure 2b rather than the composite:
  - drug names in Title Case, as the deck labels them (composite uses lower case)
  - BOLD marks compounds also significant in the FT bidirectional core
    (the composite used a star; the caption defines bold)
  - the composite's other marker, an arrow pair on vemurafenib meaning "hits both the
    INHIBIT and ACTIVATE axes", is NOT drawn: the current figure drops it and the caption
    does not define it. Set MARK_BOTH_AXES = True to restore it.
  - no panel title (the deck crops it off); the "adj. P = 0.05" line label is kept.

Writes Figure2b_crossmodel_consensus.svg (editable text, Arial) and a .png preview,
to <results>/02_geneformer_in_silico_screen/figures/.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Section path configuration (GF_WORK, results folder, shipped data): see ../gf_paths.py.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gf_paths  # noqa: E402

TSV = gf_paths.drug_table("crossmodel_consensus.tsv")
STEM = os.path.join(gf_paths.gf_results_dir("figures"), "Figure2b_crossmodel_consensus")

MARK_BOTH_AXES = False
EXTRA_ARTIFACTS = {"vx", "sodium"}      # truncated non-therapeutic perturbagen names
TITLECASE = {"vemurafenib": "Vemurafenib"}   # deck drops the "(PLX4032)" suffix

# colours and ink exactly as the original script
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
FT_C, MT_C = "#2a78d6", "#eb6834"

plt.rcParams.update({
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})

d = pd.read_csv(TSV, sep="\t")
d = d[~d["drug"].isin(EXTRA_ARTIFACTS)].copy()
d["ft"] = -np.log10(d["best_adjp_FT"])
d["mt"] = -np.log10(d["best_adjp_MT"])
d = d.sort_values("combined_signal").reset_index(drop=True)   # same order as published
y = np.arange(len(d))

fig, ax = plt.subplots(figsize=(5.2, 3.9))
ax.axvline(-np.log10(0.05), color=MUTED, lw=0.8, ls=(0, (3, 3)), zorder=0)
ax.text(-np.log10(0.05), len(d) - 0.42, " adj. P = 0.05", fontsize=7,
        color=MUTED, va="bottom")
for i, r in d.iterrows():
    ax.plot([r.ft, r.mt], [i, i], color=GRID, lw=2, zorder=1, solid_capstyle="round")
ax.scatter(d["ft"], y, s=60, color=FT_C, ec="white", lw=1.5, zorder=3,
           label="FT (excit-specific)")
ax.scatter(d["mt"], y, s=60, color=MT_C, ec="white", lw=1.5, zorder=3,
           label="MT (JointPG multi-task)")

labels = []
for _, r in d.iterrows():
    s = TITLECASE.get(r.drug, r.drug[:1].upper() + r.drug[1:])
    if MARK_BOTH_AXES and r.both_axes_either:
        s += "  ⇅"
    labels.append(s)
ax.set_yticks(y)
ax.set_yticklabels(labels)
for t, (_, r) in zip(ax.get_yticklabels(), d.iterrows()):
    if r.in_FT_core:                      # caption: bold = also in FT bidirectional core
        t.set_weight("bold")

ax.set_xlabel("−log10 adjusted-p (drug signature vs targets)")
ax.set_xlim(0, max(d.ft.max(), d.mt.max()) * 1.12)
ax.set_ylim(-0.6, len(d) - 0.1)
ax.grid(axis="x", color=GRID, lw=0.6)
ax.set_axisbelow(True)
ax.legend(handles=[Line2D([], [], marker="o", ls="none", mfc=FT_C, mec="white", mew=1.2,
                          ms=7, label="FT (excit-specific)"),
                   Line2D([], [], marker="o", ls="none", mfc=MT_C, mec="white", mew=1.2,
                          ms=7, label="MT (JointPG multi-task)")],
          loc="lower right", frameon=False, fontsize=8)

for ext in ("svg", "png"):
    fig.savefig(STEM + "." + ext, **({"dpi": 300} if ext == "png" else {}))
    print("wrote", STEM + "." + ext)

print("\nvalues plotted (from %s):" % os.path.basename(TSV))
print("  %-18s %8s %8s  %s" % ("drug", "FT", "MT", "bold (FT core)"))
for _, r in d[::-1].iterrows():
    print("  %-18s %8.3f %8.3f  %s"
          % (r.drug, r.ft, r.mt, "yes" if r.in_FT_core else ""))
