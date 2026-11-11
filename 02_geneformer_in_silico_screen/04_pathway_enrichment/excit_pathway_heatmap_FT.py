#!/usr/bin/env python
"""FT (excit-specific fine-tuned) pathway-enrichment heatmap across the FULL 2x2
perturbation matrix (both directions x both perturbation types):

    RESCUE     (TBI -> control):  delete | overexpress
    INDUCTION  (control -> TBI):  delete | overexpress

    rows  = union of each arm's TOP_PER_ARM strongest significant pathways
            (adjP<0.05), grouped by source library
    color = -log10(adjP)   (white = not enriched in that arm)

Row cap is needed because the induction-overexpress arm alone yields ~336
significant terms (its hit set is ~2.6k of a 7.4k background -> broad/low-
specificity); plotting all would be unreadable. Foreground per arm uses the
same average-shift hit definition as the MT/FT figures
(excit_pathway_heatmap_enrich.py): FDR-sig, det>=20%, shift>significant-set avg.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Section path configuration (GF_WORK, results folder, shipped data): see ../gf_paths.py.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gf_paths  # noqa: E402

# Per-arm Enrichr tables written by excit_pathway_heatmap_enrich.py: a rerun below
# <results>/02_geneformer_in_silico_screen/pathway_heatmap/enrich, else the original
# project tree $GF_WORK/results/pathway_heatmap/enrich. An arm found in neither is
# read from the shipped source matrix of this figure (see load()).
OUTDIR = gf_paths.gf_results_dir("pathway_heatmap")
OUT = os.path.join(OUTDIR, "excit_pathway_enrichment_heatmap_FT.png")
TOP_PER_ARM = 15

ARMS = [
    ("FT_delete",  "delete\nTBI→ctrl\n(rescue)"),
    ("FT_OE",      "overexpr\nTBI→ctrl\n(rescue)"),
    ("FT_del_c2t", "delete\nctrl→TBI\n(induct)"),
    ("FT_oe_c2t",  "overexpr\nctrl→TBI\n(induct)"),
]
KEYS = [a for a, _ in ARMS]
LIB_ORDER = ["GO_Biological_Process_2026", "KEGG_2026",
             "Reactome_Pathways_2024", "MSigDB_Hallmark_2020"]
LIB_SHORT = {"GO_Biological_Process_2026": "GO-BP", "KEGG_2026": "KEGG",
             "Reactome_Pathways_2024": "Reactome", "MSigDB_Hallmark_2020": "Hallmark"}
SIG = 0.05


def load(label):
    p = gf_paths.analysis_input("pathway_heatmap", "enrich", label,
                                f"{label}_enrichr_allterms.csv")
    if not os.path.exists(p):
        # No per-arm Enrichr table: take the arm's adjusted P-values from the
        # shipped matrix (data/.../pathway_enrichment/), which this script wrote
        # from those tables in the published run.
        shipped = gf_paths.shipped_pathway_arm(label)
        if shipped is not None:
            print(f"[{label}] no per-arm Enrichr table; using the shipped matrix "
                  f"{os.path.basename(gf_paths.SHIPPED_PATHWAY_MATRIX[label.split('_', 1)[0]])}")
            return shipped.rename(columns={"Adjusted P-value": label})
        print(f"[{label}] allterms MISSING ({p})")
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
    raise SystemExit("no FT enrichment results found")
mat = pd.concat(frames, axis=1)
for k in KEYS:
    if k not in mat.columns:
        mat[k] = np.nan
mat = mat[KEYS]

# rows = union of each ready arm's top-N significant terms
keep_idx = set()
for k in ready:
    s = mat[k][mat[k] < SIG].sort_values()
    keep_idx.update(s.head(TOP_PER_ARM).index.tolist())
mat = mat.loc[sorted(keep_idx)].copy()
nsig_total = {k: int((mat[k] < SIG).sum()) for k in ready}
print(f"FT arms ready: {ready}")
print(f"displayed terms (union of top{TOP_PER_ARM} sig per arm): {len(mat)}")
print(f"  sig terms shown per arm: {nsig_total}")

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
ax.set_xlabel("FT excit-specific model — perturbation arm")
ax.set_title("FT excit model, excitatory neurons — pathway enrichment across "
             "directions\n(-log10 adjP; * = adjP<0.05; left 2 = rescue TBI→ctrl, "
             f"right 2 = induction ctrl→TBI; rows = top{TOP_PER_ARM}/arm)",
             fontsize=10)

for i in range(n):
    for j, k in enumerate(KEYS):
        v = mat[k].iloc[i]
        if pd.isna(v):
            continue
        if v < SIG:
            ax.text(j, i, "*", ha="center", va="center", color="white",
                    fontsize=10, fontweight="bold")

ax.axvline(1.5, color="red", lw=1.4)   # rescue | induction divider
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

# full (uncapped) matrix dump
full = pd.concat([load(k).set_index(["Gene_set", "Term"]) for k in ready], axis=1)
full.reset_index().to_csv(
    os.path.join(OUTDIR, "excit_pathway_enrichment_matrix_FT.tsv"),
    sep="\t", index=False)
print("wrote full matrix tsv")
