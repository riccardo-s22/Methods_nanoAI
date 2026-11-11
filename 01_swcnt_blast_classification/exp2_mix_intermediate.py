"""Do the dilution mixes show an INTERMEDIATE spectral phenotype?
Avoid the whitened LDA axis (can amplify noisy directions). Use:
 (1) simple un-whitened sham->TBI centroid projection in shape space (0=sham,1=TBI);
 (2) cosine similarity of each group's (mean - sham) difference-spectrum to the TBI
     difference-spectrum  -> is it the SAME spectral change?
 (3) projection magnitude onto the TBI difference direction -> how MUCH of it (fraction).
 (4) overlay difference-from-sham spectra per group.
"""
import os
os.environ.setdefault("MIX_INVERTED", "1")  # col10=1/4 TBI, col12=3/4 TBI
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
TAG = os.environ.get("EXP2_TAG", "2h")
import numpy as np, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from repo_paths import results_dir
from exp2_load import load_all, labels, EXC
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def pr(*a): print(*a); sys.stdout.flush()

cube, em = load_all(); rows, cols, grp = labels()
def shp(C):
    S = C.copy()
    for k in range(C.shape[2]): S[:, :, k] /= S[:, :, k].sum(1, keepdims=True)
    return S.reshape(len(C), -1)

groups = ["negctrl","sham","mix25","mix50","mix75","TBI"]
frac = {"negctrl":np.nan,"sham":0.0,"mix25":0.25,"mix50":0.5,"mix75":0.75,"TBI":1.0}
def wells(g):
    idx = np.where(grp == g)[0]
    if g == "negctrl": idx = idx[idx != 0]   # drop A1 ramp
    return idx
X = {g: shp(cube[wells(g)]) for g in groups}

sham_mu = X["sham"].mean(0)
tbi_mu  = X["TBI"].mean(0)
d = tbi_mu - sham_mu                      # sham->TBI direction (un-whitened)
dn = d / (d @ d)
# (1)+(3) per-well projection: 0 at sham centroid, 1 at TBI centroid
pr("=== (1) un-whitened sham->TBI projection  (0=sham, 1=TBI) ===")
proj = {}
for g in groups:
    p = (X[g] - sham_mu) @ dn
    proj[g] = p
    pr("  %-8s frac=%4s  n=%2d  proj = %+.2f ± %.2f" % (g, frac[g], len(p), p.mean(), p.std()))
from scipy.stats import spearmanr
fr = np.array([frac[g] for g in ["sham","mix25","mix50","mix75","TBI"]])
mn = np.array([proj[g].mean() for g in ["sham","mix25","mix50","mix75","TBI"]])
pr("  monotonicity vs TBI-fraction: Spearman r=%.3f" % spearmanr(fr, mn)[0])

# (2) is each group's spectral CHANGE the same SHAPE as TBI's? cosine to TBI diff
pr("\n=== (2) cosine similarity of (group_mean - sham) to (TBI - sham) ===")
for g in ["mix25","mix50","mix75","TBI","negctrl"]:
    dg = X[g].mean(0) - sham_mu
    cos = (dg @ d) / (np.linalg.norm(dg) * np.linalg.norm(d) + 1e-30)
    pr("  %-8s cos=%+.3f   |change|/|TBI change| = %.2f" %
       (g, cos, np.linalg.norm(dg) / np.linalg.norm(d)))

# (4) overlay difference-from-sham spectra (averaged over excitations) per group
ds_by_g = {}
for g in groups:
    Sg = cube[wells(g)].copy()
    for k in range(5): Sg[:, :, k] /= Sg[:, :, k].sum(1, keepdims=True)
    Ss = cube[wells("sham")].copy()
    for k in range(5): Ss[:, :, k] /= Ss[:, :, k].sum(1, keepdims=True)
    ds_by_g[g] = (Sg.mean(0) - Ss.mean(0)).mean(1)   # (512,)

fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
colmap = {"negctrl":"#7f8c8d","sham":"#2e86c1","mix25":"#a9dfbf","mix50":"#52be80","mix75":"#1e8449","TBI":"#cb4335"}
for g in groups:
    if g == "sham": continue
    ax[0].plot(em, ds_by_g[g], label="%s (%.0f%% TBI)" % (g, 100*frac[g]) if not np.isnan(frac[g]) else g,
               c=colmap[g], lw=2 if g in ("TBI","negctrl") else 1.6)
ax[0].axhline(0, c="grey", ls="--"); ax[0].set_xlim(900,1450)
ax[0].set_xlabel("emission (nm)"); ax[0].set_ylabel("group mean − sham (shape)")
ax[0].set_title("Difference-from-sham spectra by group\n(do mixes nest between sham & TBI?)"); ax[0].legend(fontsize=8)
# right: projection vs TBI fraction, per-well
for g in ["sham","mix25","mix50","mix75","TBI"]:
    ax[1].scatter([frac[g]]*len(proj[g]), proj[g], c=colmap[g], s=26, edgecolor="k", lw=.3, zorder=3)
ax[1].plot(fr, mn, "-o", c="k", lw=1.8, zorder=2, label="group mean")
ax[1].plot([0,1],[0,1], ls=":", c="red", label="ideal linear dose")
ax[1].axhline(proj["negctrl"].mean(), c="#7f8c8d", ls="--", label="negctrl")
ax[1].set_xlabel("TBI-media fraction"); ax[1].set_ylabel("un-whitened sham→TBI projection")
ax[1].set_title("Intermediate phenotype check"); ax[1].legend(fontsize=8)
fig.tight_layout()
OUT = os.environ.get("EXP2_OUT") or results_dir("01_swcnt_blast_classification")
os.makedirs(OUT, exist_ok=True)
p = os.path.join(OUT, "exp2_%s_mix_intermediate.png" % TAG); fig.savefig(p, dpi=150)
pr("\nsaved", p)
