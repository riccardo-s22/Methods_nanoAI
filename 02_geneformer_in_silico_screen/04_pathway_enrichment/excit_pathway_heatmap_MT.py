#!/usr/bin/env python
"""MT-only pathway-enrichment heatmap for excitatory neurons across the full
2x2 perturbation matrix (both directions x both perturbation types):

    RESCUE     (TBI -> control):  delete | overexpress
    INDUCTION  (control -> TBI):  delete | overexpress

    rows  = pathways significantly enriched (adjP<0.05) in >=1 arm
    color = -log10(adjP)   (white = not enriched / arm not ready)

mean_shift>0 always means "toward the goal centroid", so a rescue-arm hit moves
TBI cells toward control while an induction-arm hit moves control cells toward
TBI. A pathway that lights up in BOTH a rescue and an induction arm is a
bidirectional driver of the excit TBI state. All arms use the identical
20%-detection / top-100 / shared-background Enrichr methodology
(excit_pathway_heatmap_enrich.py).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TBI = "/work/pi_cfallini_uri_edu/Geneformer_TBI"
ENR = f"{TBI}/results/pathway_heatmap/enrich"
OUT = f"{TBI}/results/pathway_heatmap/excit_pathway_enrichment_heatmap_MT.png"

# (arm key, 3-line column label) ordered: rescue pair then induction pair
ARMS = [
    ("MT_delete",  "delete\nTBI→ctrl\n(rescue)"),
    ("MT_OE",      "overexpr\nTBI→ctrl\n(rescue)"),
    ("MT_del_c2t", "delete\nctrl→TBI\n(induct)"),
    ("MT_oe_c2t",  "overexpr\nctrl→TBI\n(induct)"),
]
KEYS = [a for a, _ in ARMS]
LIB_ORDER = ["GO_Biological_Process_2026", "KEGG_2026",
             "Reactome_Pathways_2024", "MSigDB_Hallmark_2020"]
LIB_SHORT = {"GO_Biological_Process_2026": "GO-BP", "KEGG_2026": "KEGG",
             "Reactome_Pathways_2024": "Reactome", "MSigDB_Hallmark_2020": "Hallmark"}
SIG = 0.05


def load(label):
    p = f"{ENR}/{label}/{label}_enrichr_allterms.csv"
    if not os.path.exists(p):
        print(f"[{label}] allterms MISSING -> arm not ready ({p})")
        return None
    df = pd.read_csv(p)[["Gene_set", "Term", "Adjusted P-value"]]
    return df.rename(columns={"Adjusted P-value": label})


frames, ready = [], []
for k in KEYS:
    d = load(k)
    if d is not None:
        frames.append(d.set_index(["Gene_set", "Term"]))
        ready.append(k)
if not frames:
    raise SystemExit("no MT enrichment results found")
mat = pd.concat(frames, axis=1)
for k in KEYS:
    if k not in mat.columns:
        mat[k] = np.nan
mat = mat[KEYS]

sig_any = (mat[ready] < SIG).any(axis=1)
mat = mat[sig_any].copy()
print(f"MT arms ready: {ready}")
print(f"significant terms (adjP<{SIG} in >=1 ready arm): {len(mat)}")

mat = mat.reset_index()
mat["lib_rank"] = mat["Gene_set"].map({l: i for i, l in enumerate(LIB_ORDER)})
mat["best"] = mat[KEYS].min(axis=1)
mat = mat.sort_values(["lib_rank", "best"]).reset_index(drop=True)

L = -np.log10(mat[KEYS].astype(float).values)
masked = np.ma.masked_invalid(L)
vmax = np.nanpercentile(L[np.isfinite(L)], 98) if np.isfinite(L).any() else 3.0
vmax = max(vmax, -np.log10(SIG) + 0.5)

row_labels = [f"[{LIB_SHORT[g]:8s}] {t[:54]}" for g, t in zip(mat["Gene_set"], mat["Term"])]
n = len(mat)
fig, ax = plt.subplots(figsize=(8.5, max(4, 0.34 * n + 2)))
cmap = plt.cm.viridis.copy()
cmap.set_bad("white")
im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=0, vmax=vmax)

ax.set_xticks(range(len(ARMS)))
ax.set_xticklabels([lbl for _, lbl in ARMS], fontsize=9)
ax.set_yticks(range(n))
ax.set_yticklabels(row_labels, fontsize=7, family="monospace")
ax.set_xlabel("MT model — perturbation arm")
title = ("MT model, excitatory neurons — pathway enrichment across directions\n"
         "(-log10 adjP; * = adjP<0.05; left 2 = rescue TBI→ctrl, "
         "right 2 = induction ctrl→TBI)")
# MT_del_c2t is a genuine null (only 8 hits > avg-shift bar -> below the 20-gene
# floor; deletion rarely induces the TBI state), not a pending job.
pend = [k for k in KEYS if k not in ready]
if pend:
    title += (f"\n[blank arm: {', '.join(pend)} — too few drivers to enrich "
              "(<20 hits); biological null, not pending]")
ax.set_title(title, fontsize=10)

for i in range(n):
    for j, k in enumerate(KEYS):
        v = mat[k].iloc[i]
        if pd.isna(v):
            continue
        if v < SIG:
            ax.text(j, i, "*", ha="center", va="center", color="white",
                    fontsize=11, fontweight="bold")

# separator between rescue (cols 0-1) and induction (cols 2-3)
ax.axvline(1.5, color="red", lw=1.4)
prev = None
for i, g in enumerate(mat["Gene_set"]):
    if prev is not None and g != prev:
        ax.axhline(i - 0.5, color="black", lw=0.8)
    prev = g

cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
cb.set_label("-log10(adjusted P)")
fig.tight_layout()
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print("wrote", OUT)

mat[["Gene_set", "Term"] + KEYS].to_csv(
    f"{TBI}/results/pathway_heatmap/excit_pathway_enrichment_matrix_MT.tsv",
    sep="\t", index=False)
print("wrote matrix tsv")
