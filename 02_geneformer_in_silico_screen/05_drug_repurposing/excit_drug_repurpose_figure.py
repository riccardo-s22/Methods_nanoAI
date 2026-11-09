#!/usr/bin/env python
"""Figure: connectivity-map drug repurposing for the excitatory-neuron TBI
state (built on the FT 2x2 ISP screen; cross-checked with the MT/JointPG model).

  A  workflow schematic: 4 ISP arms -> direction-aware INHIBIT / ACTIVATE
     target sets -> Enrichr vs drug DOWN / UP signatures -> per-model ranking
     -> cross-model consensus
  B  cross-model consensus drugs: best adjusted P in FT vs MT (dumbbell)
  C  FT bidirectional-core drugs: best adjusted P, overlap size, direction

Reads only finished outputs in results/drug_repurpose/ (no Enrichr calls).
Target-set sizes in panel A come from TARGET_COUNTS below (rebuilt from the
FT ranked tables with excit detection >=20%; re-run excit_drug_repurpose.py on
Unity to print the canonical numbers and update them if they differ).

Output: results/drug_repurpose/excit_drug_repurpose_figure.{png,pdf}
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

BASE = Path(__file__).resolve().parents[1]
DR = BASE / "results/drug_repurpose"
OUT = DR / "excit_drug_repurpose_figure"

# non-therapeutic perturbagens whose truncated library names slipped past the
# ARTIFACTS list in excit_drug_repurpose.py ("VX 39793" = VX nerve agent,
# "SODIUM ARSENITE")
EXTRA_ARTIFACTS = {"vx", "sodium"}
DISPLAY = {"vemurafenib": "vemurafenib (PLX4032)"}

TARGET_COUNTS = {  # FT model, foreground genes per ISP arm and derived sets
    "ind_oe": 1429, "res_del": 988, "ind_del": 1258, "res_oe": 823,
    "inhibit": 1809, "activate": 1471, "inh_core": 594, "act_core": 596,
    "conflict": 14, "background": 7374,
}

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
FT_C, MT_C = "#2a78d6", "#eb6834"          # categorical slots 1, 2
INH_FILL, ACT_FILL, NEU_FILL = "#e8f0fb", "#fdeee6", "#f3f2ef"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
})


def box(ax, x, y, w, h, text, fill, edge=MUTED, size=8.2, weight="normal"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008,rounding_size=0.015",
                                fc=fill, ec=edge, lw=0.9))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size,
            color=INK, weight=weight, linespacing=1.35)


def arrow(ax, x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.0,
                                shrinkA=0, shrinkB=0, mutation_scale=9))


def panel_a(ax):
    c = TARGET_COUNTS
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    # column 1: ISP arms
    arms = [
        ("INDUCER · overexpress\ncontrol→TBI", c["ind_oe"], INH_FILL, 0.80),
        ("RESCUE · delete\nTBI→control", c["res_del"], INH_FILL, 0.585),
        ("INDUCER · delete\ncontrol→TBI", c["ind_del"], ACT_FILL, 0.315),
        ("RESCUE · overexpress\nTBI→control", c["res_oe"], ACT_FILL, 0.10),
    ]
    ax.text(0.095, 0.985, "Fine-tuned excit ISP arms", ha="center", va="top",
            fontsize=8.5, color=INK2, weight="bold")
    for lab, n, fill, y in arms:
        box(ax, 0.012, y, 0.175, 0.13, f"{lab}\n{n:,} driver genes", fill)
    # column 2: target sets
    ax.text(0.35, 0.985, "Direction-aware targets", ha="center", va="top",
            fontsize=8.5, color=INK2, weight="bold")
    box(ax, 0.255, 0.60, 0.19, 0.25,
        f"INHIBIT targets\n(raising gene → TBI;\nlowering it rescues)\n\n"
        f"union {c['inhibit']:,}\nbidirectional core {c['inh_core']}", INH_FILL)
    box(ax, 0.255, 0.115, 0.19, 0.25,
        f"ACTIVATE targets\n(losing gene → TBI;\nraising it rescues)\n\n"
        f"union {c['activate']:,}\nbidirectional core {c['act_core']}", ACT_FILL)
    for _, _, _, y in arms[:2]:
        arrow(ax, 0.187, y + 0.065, 0.253, 0.725)
    for _, _, _, y in arms[2:]:
        arrow(ax, 0.187, y + 0.065, 0.253, 0.24)
    ax.text(0.35, 0.475, f"{c['conflict']} direction-conflicting\ngenes dropped",
            ha="center", va="center", fontsize=7.5, color=MUTED, style="italic")
    # column 3: connectivity query
    ax.text(0.585, 0.985, "Connectivity query (Enrichr)", ha="center", va="top",
            fontsize=8.5, color=INK2, weight="bold")
    box(ax, 0.49, 0.62, 0.19, 0.21,
        "drug DOWN signatures\nLINCS L1000 + GEO\n\ndrug must LOWER\nINHIBIT targets", NEU_FILL)
    box(ax, 0.49, 0.135, 0.19, 0.21,
        "drug UP signatures\nLINCS L1000 + GEO\n\ndrug must RAISE\nACTIVATE targets", NEU_FILL)
    arrow(ax, 0.447, 0.725, 0.488, 0.725)
    arrow(ax, 0.447, 0.24, 0.488, 0.24)
    ax.text(0.585, 0.475, f"background: {c['background']:,} genes\n"
            "detected in ≥20% of excit cells\nBH adj. P < 0.05",
            ha="center", va="center", fontsize=7.5, color=MUTED, style="italic")
    # column 4: ranking + consensus
    ax.text(0.85, 0.985, "Ranking & consensus", ha="center", va="top",
            fontsize=8.5, color=INK2, weight="bold")
    box(ax, 0.73, 0.555, 0.25, 0.33,
        "per-drug score = Σ −log10 adj. P\n(toxicant / chemo signatures removed)\n\n"
        "FT model: 45 drugs (union)\n6 drugs (bidirectional core)\n"
        "MT (JointPG) model: 21 drugs", NEU_FILL)
    box(ax, 0.73, 0.135, 0.25, 0.25,
        "Cross-model consensus\n(significant in FT AND MT)\n\n"
        "6 drugs  ·  2 also in FT core:\nresveratrol, diclofenac",
        "#ffffff", edge=INK, weight="bold")
    arrow(ax, 0.682, 0.725, 0.728, 0.725)
    arrow(ax, 0.682, 0.24, 0.728, 0.66)
    arrow(ax, 0.855, 0.553, 0.855, 0.387)


def panel_b(ax):
    d = pd.read_csv(DR / "crossmodel_consensus.tsv", sep="\t")
    d = d[~d["drug"].isin(EXTRA_ARTIFACTS)].copy()
    d["ft"] = -np.log10(d["best_adjp_FT"])
    d["mt"] = -np.log10(d["best_adjp_MT"])
    d = d.sort_values("combined_signal").reset_index(drop=True)
    y = np.arange(len(d))
    ax.axvline(-np.log10(0.05), color=MUTED, lw=0.8, ls=(0, (3, 3)), zorder=0)
    ax.text(-np.log10(0.05), len(d) - 0.35, " adj. P = 0.05", fontsize=7,
            color=MUTED, va="bottom")
    for i, r in d.iterrows():
        ax.plot([r.ft, r.mt], [i, i], color=GRID, lw=2, zorder=1,
                solid_capstyle="round")
    ax.scatter(d["ft"], y, s=60, color=FT_C, ec="white", lw=1.5, zorder=3,
               label="FT (excit-specific)")
    ax.scatter(d["mt"], y, s=60, color=MT_C, ec="white", lw=1.5, zorder=3,
               label="MT (JointPG multi-task)")
    labels = []
    for _, r in d.iterrows():
        s = DISPLAY.get(r.drug, r.drug)
        if r.in_FT_core:
            s += "  ★"
        if r.both_axes_either:
            s += "  ⇅"
        labels.append(s)
    ax.set_yticks(y, labels)
    for t, (_, r) in zip(ax.get_yticklabels(), d.iterrows()):
        if r.in_FT_core:
            t.set_weight("bold")
    ax.set_xlabel("best −log10 adjusted P (drug signature vs targets)")
    ax.set_xlim(0, max(d.ft.max(), d.mt.max()) * 1.12)
    ax.set_ylim(-0.6, len(d) - 0.1)
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.legend(loc="lower right", frameon=False, fontsize=8)
    ax.set_title("Cross-model consensus drugs\n★ also in FT bidirectional core   "
                 "⇅ hits both INHIBIT and ACTIVATE axes", loc="left",
                 fontsize=9, color=INK)


def panel_c(ax):
    d = pd.read_csv(DR / "FT/drug_candidates_ranked_core.tsv", sep="\t")
    d = d[~d["artifact"]].copy()
    t = pd.read_csv(DR / "FT/all_sig_drug_terms.tsv", sep="\t")
    t["n_genes"] = t["Genes"].str.split(";").str.len()
    t["src"] = t["Term"].str.extract(r"(GSE\d+)")[0] + " · " + \
        t["Term"].str.extract(r"(human|mouse|rat)")[0]
    best = t.sort_values("Adjusted P-value").drop_duplicates("drug")
    d = d.merge(best[["drug", "n_genes", "src"]], on="drug", how="left")
    d["score"] = -np.log10(d["best_adjp"])
    d = d.sort_values("score").reset_index(drop=True)
    y = np.arange(len(d))
    ax.barh(y, d["score"], height=0.55, color=FT_C, ec="white", lw=2)
    ax.axvline(-np.log10(0.05), color=MUTED, lw=0.8, ls=(0, (3, 3)), zorder=0)
    for i, r in d.iterrows():
        arrow_txt = "↓ lowers INHIBIT targets" if r.directions.startswith("down") \
            else "↑ raises ACTIVATE targets"
        ax.text(r.score + 0.12, i, f"{arrow_txt}  ·  {int(r.n_genes)} genes  ·  {r.src}",
                va="center", fontsize=7.3, color=INK2)
    ax.set_yticks(y, d["drug"])
    ax.set_xlim(0, d["score"].max() * 2.3)
    ax.set_xlabel("−log10 adjusted P")
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.set_title("FT bidirectional-core drugs (targets confirmed\nby BOTH induction and rescue screens)",
                 loc="left", fontsize=9, color=INK)


def main():
    fig = plt.figure(figsize=(15, 10.2), facecolor="white")
    gs = fig.add_gridspec(2, 2, height_ratios=[0.8, 1.0], width_ratios=[1, 1.05],
                          hspace=0.22, wspace=0.42)
    axa = fig.add_subplot(gs[0, :]); panel_a(axa)
    axb = fig.add_subplot(gs[1, 0]); panel_b(axb)
    axc = fig.add_subplot(gs[1, 1]); panel_c(axc)
    for ax, letter, dx in [(axa, "A", -0.01), (axb, "B", -0.27), (axc, "C", -0.22)]:
        ax.text(dx, 1.06 if letter == "A" else 1.16, letter, transform=ax.transAxes,
                fontsize=15, weight="bold", color=INK, va="top")
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}.{ext}", dpi=300, bbox_inches="tight")
        print(f"wrote {OUT}.{ext}")


if __name__ == "__main__":
    main()
