"""exp2/2h: (1) project negctrl + dilution mixes onto sham-vs-TBI axis (dose-response,
   independent of even/odd parity since mixes are single columns 10/11/12);
   (2) identify emission peaks driving the shape difference."""
import os
os.environ.setdefault("MIX_INVERTED", "1")  # col10=1/4 TBI, col12=3/4 TBI
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
TAG = os.environ.get("EXP2_TAG", "2h")
import numpy as np, sys, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from repo_paths import results_dir
from exp2_load import load_all, labels, EXC
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from scipy.stats import mannwhitneyu
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def pr(*a): print(*a); sys.stdout.flush()

cube, em = load_all(); rows, cols, grp = labels()
valid = grp != "negctrl"  # keep all; but A1 is negctrl & ramp -> handled (index0 grp negctrl)
def shape(C):
    S = C.copy()
    for k in range(C.shape[2]): S[:, :, k] /= S[:, :, k].sum(1, keepdims=True)
    return S.reshape(len(C), -1)

mask = np.isin(grp, ["sham", "TBI"])
y = (grp[mask] == "TBI").astype(int)
Xtr = shape(cube[mask])
pipe = make_pipeline(StandardScaler(), PCA(20, random_state=0)).fit(Xtr)
Ztr = pipe.transform(Xtr)
lda = LDA().fit(Ztr, y)

def score(groupname):
    m = grp == groupname
    if groupname == "negctrl":  # drop A1 ramp (index0)
        idx = np.where(m)[0]; idx = idx[idx != 0]; m = np.zeros(96, bool); m[idx] = True
    Z = pipe.transform(shape(cube[m]))
    return lda.decision_function(Z)

pr("=== sham-vs-TBI LDA decision score by group (higher = more TBI-like) ===")
order = ["negctrl", "sham", "mix25", "mix50", "mix75", "TBI"]
frac = {"negctrl": np.nan, "sham": 0.0, "mix25": 0.25, "mix50": 0.5, "mix75": 0.75, "TBI": 1.0}
means = {}
for g in order:
    s = score(g); means[g] = s.mean()
    pr("  %-8s (TBI frac %4s)  n=%2d  score %+.2f ± %.2f" %
       (g, frac[g], len(s), s.mean(), s.std()))
# dose monotonicity across mixes (exclude negctrl)
fr = np.array([frac[g] for g in ["sham","mix25","mix50","mix75","TBI"]])
mn = np.array([means[g] for g in ["sham","mix25","mix50","mix75","TBI"]])
from scipy.stats import spearmanr, kendalltau
pr("\n  dose monotonicity (sham..TBI): Spearman r=%.3f  Kendall tau=%.3f" %
   (spearmanr(fr, mn)[0], kendalltau(fr, mn)[0]))
# do the intermediate mixes differ from sham? from TBI?
shamS, tbiS = score("sham"), score("TBI")
for g in ["mix25","mix50","mix75"]:
    s = score(g)
    pr("  %s vs sham MWU p=%.3f | vs TBI MWU p=%.3f" %
       (g, mannwhitneyu(s, shamS).pvalue, mannwhitneyu(s, tbiS).pvalue))

# ---- (2) driver emission peaks: per-emission |TBI-sham| of shape spectra, averaged over excitations
pr("\n=== driver emission regions (shape difference magnitude) ===")
Csham, Ctbi = cube[grp=="sham"], cube[grp=="TBI"]
def shp(C):
    S=C.copy()
    for k in range(C.shape[2]): S[:,:,k]/=S[:,:,k].sum(1,keepdims=True)
    return S
ds = shp(Ctbi).mean(0) - shp(Csham).mean(0)   # (512, 5)
dmag = np.abs(ds).mean(1)
peaks_em = em[np.argsort(dmag)[::-1][:8]]
pr("  top-8 emission wavelengths by |TBI-sham| shape diff: " +
   ", ".join("%.0f" % w for w in np.sort(peaks_em)))

# figure: dose-response + driver spectrum
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
xs = [frac[g] for g in ["sham","mix25","mix50","mix75","TBI"]]
ys = [means[g] for g in ["sham","mix25","mix50","mix75","TBI"]]
for g in ["sham","mix25","mix50","mix75","TBI"]:
    s = score(g); ax[0].scatter([frac[g]]*len(s), s, alpha=.5, c="#34495e", s=20)
ax[0].plot(xs, ys, "-o", c="#cb4335", lw=2, label="group mean")
ax[0].axhline(score("negctrl").mean(), c="grey", ls=":", label="negctrl")
ax[0].set_xlabel("TBI-media fraction"); ax[0].set_ylabel("sham-vs-TBI LDA score")
ax[0].set_title("Dilution dose-response (trained on sham/TBI only)"); ax[0].legend()
for k, ex in enumerate(EXC):
    ax[1].plot(em, ds[:, k], lw=1.2, label="ex%d" % ex)
ax[1].plot(em, ds.mean(1), c="k", lw=2.2, label="mean")
ax[1].axhline(0, c="grey", ls="--"); ax[1].set_xlim(900, 1450)
ax[1].set_xlabel("emission (nm)"); ax[1].set_ylabel("TBI - sham (shape)")
ax[1].set_title("Driver shape-difference spectrum"); ax[1].legend(fontsize=7)
fig.tight_layout()
OUT = os.environ.get("EXP2_OUT") or results_dir("01_swcnt_blast_classification")
os.makedirs(OUT, exist_ok=True)
p = os.path.join(OUT, "exp2_%s_dose_drivers.png" % TAG); fig.savefig(p, dpi=150)
pr("\nsaved", p)
