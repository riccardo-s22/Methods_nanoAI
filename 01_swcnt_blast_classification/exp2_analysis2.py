"""exp2/2h: permutation null (LDA, fast) + per-excitation + broad-dimming + peak-shift + figure."""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"): os.environ[v]="1"
TAG = os.environ.get("EXP2_TAG", "2h")
import numpy as np, warnings, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from repo_paths import results_dir
warnings.filterwarnings("ignore")
from exp2_load import load_all, labels, EXC
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut, LeaveOneGroupOut
from sklearn.metrics import roc_auc_score
from scipy.stats import mannwhitneyu
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def pr(*a): print(*a); sys.stdout.flush()

cube, em = load_all()
rows, cols, grp = labels()
mask = np.isin(grp, ["sham", "TBI"])
y = (grp[mask] == "TBI").astype(int); rr = rows[mask]; C = cube[mask]
Sn = C.copy()
for k in range(C.shape[2]): Sn[:, :, k] /= Sn[:, :, k].sum(1, keepdims=True)
Xs = Sn.reshape(len(C), -1)
Z = make_pipeline(StandardScaler(), PCA(20, random_state=0)).fit_transform(Xs)

# LEGACY (kept for comparison only): PCA fitted on ALL wells before CV -> leaks the
# held-out rows into the feature basis; AUC read from predict_proba, which saturates
# at exactly 1.0 for LDA and makes roc_auc_score charge tied pairs 0.5.
def row_out_auc_legacy(Z, yy, groups, clf):
    p = np.zeros(len(yy))
    for tr, te in LeaveOneGroupOut().split(Z, yy, groups):
        m = clf().fit(Z[tr], yy[tr])
        p[te] = m.predict_proba(Z[te])[:, 1] if hasattr(m, "predict_proba") else m.decision_function(Z[te])
    return roc_auc_score(yy, p)

# CORRECTED: scaler + PCA refitted inside every fold, score from decision_function.
def row_out_auc(X, yy, groups, clf, npca=20):
    p = np.zeros(len(yy))
    for tr, te in LeaveOneGroupOut().split(X, yy, groups):
        m = make_pipeline(StandardScaler(), PCA(npca, random_state=0), clf()).fit(X[tr], yy[tr])
        p[te] = m.decision_function(X[te])
    return roc_auc_score(yy, p)

pr("=== permutation null (row-out CV, LDA on shape spectra, n=1000) ===")
obs_leg = row_out_auc_legacy(Z, y, rr, LDA)
rng = np.random.default_rng(0); c = 0; n = 1000
for _ in range(n):
    if row_out_auc_legacy(Z, rng.permutation(y), rr, LDA) >= obs_leg: c += 1
pr("  [legacy: global PCA20 + predict_proba]  AUC=%.4f  perm p=%.4f" % (obs_leg, (c+1)/(n+1)))

obs_lda = row_out_auc(Xs, y, rr, LDA)
rng = np.random.default_rng(0); c = 0
for _ in range(n):
    if row_out_auc(Xs, rng.permutation(y), rr, LDA) >= obs_lda: c += 1
pr("  [corrected: in-fold PCA20 + decision_function]  AUC=%.4f  perm p=%.4f" % (obs_lda, (c+1)/(n+1)))

pr("\n=== per-excitation (shape-normalized, LOO PCA8 LogReg) ===")
def loo(X, yy, npca=8):
    p = np.zeros(len(yy))
    for tr, te in LeaveOneOut().split(X):
        m = make_pipeline(StandardScaler(), PCA(npca, random_state=0), LogisticRegression(max_iter=2000))
        m.fit(X[tr], yy[tr]); p[te] = m.predict_proba(X[te])[:, 1]
    return roc_auc_score(yy, p)
per_ex = []
for k, ex in enumerate(EXC):
    Xk = C[:, :, k].copy(); Xk = Xk / Xk.sum(1, keepdims=True)
    a = loo(Xk, y); per_ex.append(a); pr("  ex%3d  LOO AUC=%.3f" % (ex, a))

pr("\n=== broad-dimming + emission peak shift ===")
tot = C.sum((1, 2))
pr("  total intensity: sham %.3g TBI %.3g  ratio %.3f  MWU p=%.2g" %
   (tot[y==0].mean(), tot[y==1].mean(), tot[y==1].mean()/tot[y==0].mean(), mannwhitneyu(tot[y==0], tot[y==1]).pvalue))
cents = []
for k, ex in enumerate(EXC):
    sp = C[:, :, k]; cen = (sp * em).sum(1) / sp.sum(1); cents.append((cen[y==0].mean(), cen[y==1].mean()))
    pr("  ex%3d centroid sham %.2f TBI %.2f  shift %+.2f nm  MWU p=%.2g" %
       (ex, cen[y==0].mean(), cen[y==1].mean(), cen[y==1].mean()-cen[y==0].mean(), mannwhitneyu(cen[y==0], cen[y==1]).pvalue))

# ---- figure ----
OUT = os.environ.get("EXP2_OUT") or results_dir("01_swcnt_blast_classification")
os.makedirs(OUT, exist_ok=True)
fig, ax = plt.subplots(2, 2, figsize=(12, 8))
# (0,0) mean emission spectra at ex670 (brightest), sham vs TBI, shape-normalized
k0 = EXC.index(670)
sp = Sn[:, :, k0]
ax[0,0].plot(em, sp[y==0].mean(0), label="sham", c="#2e86c1", lw=2)
ax[0,0].plot(em, sp[y==1].mean(0), label="TBI", c="#cb4335", lw=2)
ax[0,0].set_xlim(900, 1450); ax[0,0].set_xlabel("emission (nm)"); ax[0,0].set_ylabel("norm. intensity")
ax[0,0].set_title("Mean emission spectrum @ex670 (shape-normalized)"); ax[0,0].legend()
# (0,1) difference spectrum across excitations
for k, ex in enumerate(EXC):
    d = Sn[:, :, k][y==1].mean(0) - Sn[:, :, k][y==0].mean(0)
    ax[0,1].plot(em, d, label="ex%d" % ex, lw=1.3)
ax[0,1].axhline(0, c="grey", ls="--"); ax[0,1].set_xlim(900, 1450)
ax[0,1].set_xlabel("emission (nm)"); ax[0,1].set_ylabel("TBI - sham (norm)")
ax[0,1].set_title("Difference spectra (TBI - sham)"); ax[0,1].legend(fontsize=7)
# (1,0) per-excitation AUC + total-dimming
ax[1,0].bar(range(len(EXC)), per_ex, color="#5dade2", edgecolor="k")
ax[1,0].axhline(0.5, c="grey", ls="--"); ax[1,0].set_ylim(0, 1.05)
ax[1,0].set_xticks(range(len(EXC))); ax[1,0].set_xticklabels(EXC)
ax[1,0].set_xlabel("excitation (nm)"); ax[1,0].set_ylabel("LOO AUC")
ax[1,0].set_title("Per-excitation sham-vs-TBI separability")
# (1,1) PCA scatter (first 2 LDA-relevant comps) colored by group
zz = make_pipeline(StandardScaler(), PCA(2, random_state=0)).fit_transform(Xs)
ax[1,1].scatter(zz[y==0,0], zz[y==0,1], c="#2e86c1", label="sham", s=40, edgecolor="k", lw=.3)
ax[1,1].scatter(zz[y==1,0], zz[y==1,1], c="#cb4335", label="TBI", s=40, edgecolor="k", lw=.3)
ax[1,1].set_xlabel("PC1"); ax[1,1].set_ylabel("PC2")
ax[1,1].set_title("Shape-spectrum PCA (sham vs TBI)"); ax[1,1].legend()
fig.suptitle("exp2 / %s  —  interleaved layout, sham vs TBI" % TAG, fontsize=13)
fig.tight_layout(rect=[0,0,1,0.97])
pth = os.path.join(OUT, "exp2_%s_sham_tbi.png" % TAG); fig.savefig(pth, dpi=150)
pr("\nsaved figure", pth)
