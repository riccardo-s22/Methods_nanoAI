#!/usr/bin/env python
"""Heatmap of per-cell cosine shifts (toward control) for top rescue-driver
genes in the FINE-TUNED EXCITATORY-NEURON-SPECIFIC model (Geneformer_TBI),
delete arm (TBI->control). Parallel to the MT version
(Geneformer_JointPG/scripts/heatmap_excit_cosine_shifts.py).

Rows = top-N driver genes by mean_shift (eligible det>=20% via the shared
layer3 detection table, positive shift). Columns = excitatory-neuron TBI cells.
Value = cosine shift to control when the gene is deleted in that cell.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/work/pi_cfallini_uri_edu/Geneformer_TBI/scripts")
from perturb_safeguards import (make_parse_full_run, iter_per_cell_shifts, RunSpec)

TBI = "/work/pi_cfallini_uri_edu/Geneformer_TBI"
RUN = Path(f"{TBI}/results/perturb_neurons_finetuned_tbi_to_control")
CHUNKS = RUN / "chunks_delete"
UNF = RUN / "stats_delete/excit_finetuned_rescue_delete_unfiltered_all_gene_ranked.tsv"
DET = ("/work/pi_cfallini_uri_edu/Geneformer_JointPG/results/isp/"
       "garza_disease_binary_TBI_to_control/mechano_analysis/layer3/"
       "garza_qc4_de_by_celltype.tsv")
CELLTYPE_PREFIX = "excit_finetuned_rescue"
CELLS_PER_CHUNK = 56
TOTAL_START = 2047
THRESH = 0.20
TOPN = 40
OUT = str(RUN / "stats_delete/excit_finetuned_cosine_shift_heatmap.png")

# eligibility via shared layer3 detection (excitatory_neuron)
det = pd.read_csv(DET, sep="\t")
det = det[det.cell_type == "excitatory_neuron"][["ensembl_id", "det_TBI", "det_ctrl"]].copy()
det["det_max"] = det[["det_TBI", "det_ctrl"]].max(axis=1)

unf = pd.read_csv(UNF, sep="\t")
m = unf.merge(det, on="ensembl_id", how="inner")
elig = m[(m["det_max"] >= THRESH) & (m["mean_shift"] > 0)]
top = elig.sort_values("mean_shift", ascending=False).head(TOPN)
gene_tokens = list(top["gene_token"].astype(int))
tok2sym = {int(t): s for t, s in zip(top["gene_token"], top["gene_symbol"])}
tokset = set(gene_tokens)
print(f"top{TOPN} driver genes: {', '.join(top['gene_symbol'].head(15))} ...")

# all perturbed cells are excit TBI (dataset is excit-only); build full matrix
M = np.zeros((len(gene_tokens), TOTAL_START), dtype=np.float32)
row_of = {t: i for i, t in enumerate(gene_tokens)}
seen = np.zeros(TOTAL_START, dtype=bool)

spec = RunSpec(chunks_dir=CHUNKS,
               parse_fn=make_parse_full_run(CELLTYPE_PREFIX, isp_mode="delete"),
               cells_per_chunk=CELLS_PER_CHUNK, total_tbi_cells=TOTAL_START,
               goal_state="control")
n_entries = 0
for cell_idx, token, shift in iter_per_cell_shifts(spec):
    seen[cell_idx] = True
    if token in tokset:
        M[row_of[token], cell_idx] = shift
        n_entries += 1
# keep only cells actually perturbed
M = M[:, seen]
print(f"cells perturbed: {seen.sum()} / {TOTAL_START}; "
      f"filled {n_entries:,} (gene,cell) entries; matrix {M.shape}")

try:
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist
    col_order = leaves_list(linkage(pdist(M.T), method="ward"))
except Exception as e:
    print("col clustering skipped:", e)
    col_order = np.arange(M.shape[1])
Mo = M[:, col_order]

vmax = np.percentile(np.abs(M[M != 0]), 99) if (M != 0).any() else 0.01
fig, ax = plt.subplots(figsize=(16, 11))
im = ax.imshow(Mo, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
               interpolation="nearest")
ax.set_yticks(range(len(gene_tokens)))
ax.set_yticklabels([tok2sym[t] for t in gene_tokens], fontsize=8)
ax.set_xlabel(f"excitatory-neuron TBI cells (n={M.shape[1]}, ward-clustered)")
ax.set_ylabel(f"top-{TOPN} rescue-driver genes (by mean shift)")
ax.set_title("Fine-tuned excitatory-neuron model — per-cell cosine shift toward "
             "control (delete)\n(red = deletion rescues cell toward control; "
             "white = gene not expressed in cell)")
cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
cb.set_label("cosine shift to control")
fig.tight_layout()
fig.savefig(OUT, dpi=140)
print("wrote", OUT)
np.savez_compressed(str(RUN / "stats_delete/excit_finetuned_cosine_shift_matrix.npz"),
                    matrix=M, genes=np.array([tok2sym[t] for t in gene_tokens]))
print("wrote matrix npz")
