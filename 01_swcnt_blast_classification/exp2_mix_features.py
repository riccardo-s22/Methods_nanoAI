"""Does the intermediate (dose-graded) phenotype emerge in ANY classification feature?
For every shape feature: (a) sham-vs-TBI discrimination AUC; (b) dose-monotonicity =
Spearman(TBI-fraction, feature) over sham+mixes+TBI wells, sign-aligned to TBI direction.
Also scan the 20 PCA components actually fed to the classifier.
A feature in the discriminative(>0.8) AND monotonic(>0.4) quadrant = a true dose-graded marker.
"""
import os
os.environ.setdefault("MIX_INVERTED", "1")  # col10=1/4 TBI, col12=3/4 TBI
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
TAG = os.environ.get("EXP2_TAG", "2h")
import numpy as np, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from repo_paths import results_dir
from exp2_load import load_all, labels, EXC
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def pr(*a): print(*a); sys.stdout.flush()

cube, em = load_all(); rows, cols, grp = labels()
def shp(C):
    S = C.copy()
    for k in range(C.shape[2]): S[:, :, k] /= S[:, :, k].sum(1, keepdims=True)
    return S
def wells(g):
    idx = np.where(grp == g)[0]
    if g == "negctrl": idx = idx[idx != 0]
    return idx

frac = {"sham":0.0,"mix25":0.25,"mix50":0.5,"mix75":0.75,"TBI":1.0}
order = ["sham","mix25","mix50","mix75","TBI"]
idx_all = np.concatenate([wells(g) for g in order])
fr_all  = np.concatenate([np.full(len(wells(g)), frac[g]) for g in order])
y_st    = np.concatenate([np.full(len(wells(g)), 0 if g=="sham" else (1 if g=="TBI" else -1)) for g in order])

# ---- per-emission features (averaged over the 5 excitations) ----
Sall = shp(cube)                  # (96,512,5)
F = Sall.mean(2)                  # (96,512) shape per emission
Xall = F[idx_all]                 # (88,512)
shamM = (y_st==0); tbiM = (y_st==1)
disc = np.zeros(512); mono = np.zeros(512)
for j in range(512):
    v = Xall[:, j]
    sgn = np.sign(v[tbiM].mean() - v[shamM].mean()) or 1.0
    vv = v * sgn
    disc[j] = roc_auc_score(np.r_[np.zeros(shamM.sum()), np.ones(tbiM.sum())],
                            np.r_[vv[shamM], vv[tbiM]])
    mono[j] = spearmanr(fr_all, vv)[0]
pr("=== per-emission shape features (avg over excitations) ===")
gooddisc = disc > 0.85
pr("  emission bands with sham-vs-TBI AUC>0.85: %d / 512" % gooddisc.sum())
both = gooddisc & (mono > 0.4)
pr("  of those, ALSO dose-monotonic (Spearman>0.4): %d" % both.sum())
if both.sum():
    pr("  -> emission nm: " + ", ".join("%.0f(mono%.2f)"%(em[j],mono[j]) for j in np.where(both)[0][:12]))
pr("  among AUC>0.85 bands: mean dose-Spearman = %.3f  (max %.3f, min %.3f)"
   % (mono[gooddisc].mean(), mono[gooddisc].max(), mono[gooddisc].min()))
# the single most discriminative band: its group means
jbest = np.argmax(disc)
sgn = np.sign(F[wells('TBI')][:,jbest].mean()-F[wells('sham')][:,jbest].mean()) or 1
pr("  most-discriminative band %.0f nm (AUC %.2f): group means(TBI-dir) " % (em[jbest], disc[jbest])
   + " ".join("%s=%.3g"%(g, sgn*F[wells(g)][:,jbest].mean()) for g in order))

# ---- PCA components actually used by the classifier ----
pipe = make_pipeline(StandardScaler(), PCA(20, random_state=0)).fit(Sall.reshape(96,-1)[np.concatenate([wells('sham'),wells('TBI')])])
Z = pipe.transform(Sall.reshape(96,-1))
Za = Z[idx_all]
pr("\n=== 20 PCA components (classifier inputs) ===")
pr("  %3s %8s %10s" % ("PC","shamTBI_AUC","dose_Spear"))
pc_disc=[]; pc_mono=[]
for c in range(20):
    v = Za[:, c]
    sgn = np.sign(v[tbiM].mean()-v[shamM].mean()) or 1.0
    vv=v*sgn
    a = roc_auc_score(np.r_[np.zeros(shamM.sum()),np.ones(tbiM.sum())], np.r_[vv[shamM],vv[tbiM]])
    m = spearmanr(fr_all, vv)[0]
    pc_disc.append(a); pc_mono.append(m)
    if a>0.7 or abs(m)>0.4:
        gm=" ".join("%s=%.2f"%(g, sgn*Za[fr_all==frac[g]][:,c].mean()) for g in order)
        pr("  PC%-2d  %6.2f   %+7.2f   | %s" % (c, a, m, gm))

# figure: discrimination vs dose-monotonicity for all emission features + PCs
fig, ax = plt.subplots(1,2, figsize=(13,5))
sc=ax[0].scatter(disc, mono, c=em, cmap="viridis", s=10)
ax[0].axvline(0.85, c="r", ls="--"); ax[0].axhline(0.4, c="r", ls="--")
ax[0].axhline(0, c="grey", lw=.6)
ax[0].set_xlabel("sham-vs-TBI AUC (discrimination)"); ax[0].set_ylabel("dose monotonicity (Spearman vs TBI frac)")
ax[0].set_title("Per-emission features\n(top-right red box = discriminative AND dose-graded)")
plt.colorbar(sc, ax=ax[0], label="emission nm")
ax[1].scatter(pc_disc, pc_mono, c="#34495e", s=40)
for c in range(20):
    if pc_disc[c]>0.7 or abs(pc_mono[c])>0.4: ax[1].annotate("PC%d"%c,(pc_disc[c],pc_mono[c]),fontsize=8)
ax[1].axvline(0.85,c="r",ls="--"); ax[1].axhline(0.4,c="r",ls="--"); ax[1].axhline(0,c="grey",lw=.6)
ax[1].set_xlabel("sham-vs-TBI AUC"); ax[1].set_ylabel("dose monotonicity (Spearman)")
ax[1].set_title("PCA components (classifier inputs)")
fig.tight_layout()
OUT = os.environ.get("EXP2_OUT") or results_dir("01_swcnt_blast_classification")
os.makedirs(OUT, exist_ok=True)
p=os.path.join(OUT,"exp2_%s_mix_feature_scan.png" % TAG); fig.savefig(p,dpi=150)
pr("\nsaved", p)
