"""
exp2/2h sham-vs-TBI analysis. Interleaved layout decouples condition from column
position (the randomized-layout control exp1 lacked). 5 excitations, 512 emission pts.

Tests:
 A) LOO classification (raw spectra, shape-normalized spectra) PCA+LogReg / LinearSVM
 B) leave-one-ROW-out (row = dominant brightness/batch axis)
 C) label-permutation null
 D) per-excitation separability
 E) emission peak-position shift (the exp1-2h marker) + broad-dimming
 F) interleave control: even-vs-odd among NON-organoid? n/a; instead column-parity vs row
"""
import numpy as np, warnings
warnings.filterwarnings("ignore")
from exp2_load import load_all, labels, EXC
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.model_selection import LeaveOneOut, LeaveOneGroupOut
from sklearn.metrics import roc_auc_score
from scipy.stats import mannwhitneyu

cube, em = load_all()
rows, cols, grp = labels()
mask = np.isin(grp, ["sham", "TBI"])
y = (grp[mask] == "TBI").astype(int)
rr = rows[mask]
C = cube[mask]  # (64, 512, 5)
print("n sham=%d TBI=%d" % ((y==0).sum(), (y==1).sum()))

def feats(kind):
    X = C.reshape(len(C), -1).copy()  # raw concat
    if kind == "raw":
        return X
    if kind == "shape":  # L1-normalize each excitation spectrum per well (remove brightness)
        Xn = C.copy()
        for k in range(C.shape[2]):
            s = Xn[:, :, k].sum(1, keepdims=True)
            Xn[:, :, k] = Xn[:, :, k] / s
        return Xn.reshape(len(C), -1)

def loo_auc(X, yy, npca=10, clf="lr"):
    pr = np.zeros(len(yy))
    for tr, te in LeaveOneOut().split(X):
        if clf == "lr":
            m = make_pipeline(StandardScaler(), PCA(npca, random_state=0), LogisticRegression(max_iter=5000))
            m.fit(X[tr], yy[tr]); pr[te] = m.predict_proba(X[te])[:, 1]
        else:
            m = make_pipeline(StandardScaler(), PCA(npca, random_state=0), LinearSVC(C=0.1, max_iter=5000, random_state=0))
            m.fit(X[tr], yy[tr]); pr[te] = m.decision_function(X[te])
    return roc_auc_score(yy, pr)

def row_out_auc(X, yy, groups, npca=10):
    pr = np.zeros(len(yy))
    for tr, te in LeaveOneGroupOut().split(X, yy, groups):
        m = make_pipeline(StandardScaler(), PCA(npca, random_state=0), LogisticRegression(max_iter=5000))
        m.fit(X[tr], yy[tr]); pr[te] = m.predict_proba(X[te])[:, 1]
    return roc_auc_score(yy, pr)

def perm_p(X, yy, obs, fn, n=500, seed=0, **kw):
    rng = np.random.default_rng(seed); c = 0
    for _ in range(n):
        if fn(X, rng.permutation(yy), **kw) >= obs: c += 1
    return (c+1)/(n+1)

import sys
def pr(*a): print(*a); sys.stdout.flush()

pr("\n=== A) LOO classification ===")
res = {}
for kind in ["raw", "shape"]:
    X = feats(kind)
    a_lr = loo_auc(X, y, 10, "lr")
    a_sv = loo_auc(X, y, 10, "svm")
    res[kind] = (X, a_lr)
    pr("  %-6s PCA10 LogReg AUC=%.3f | LinSVM AUC=%.3f" % (kind, a_lr, a_sv))

pr("\n=== B) leave-one-ROW-out (generalize across plate rows) ===")
for kind in ["raw", "shape"]:
    X = res[kind][0]
    a = row_out_auc(X, y, rr, 10)
    pr("  %-6s row-out AUC=%.3f" % (kind, a))

pr("\n=== C) label-permutation null (n=100, shape PCA10 LOO) ===")
Xs, a_obs = res["shape"]
p = perm_p(Xs, y, a_obs, loo_auc, n=100, npca=10, clf="lr")
pr("  shape LOO AUC=%.3f  perm p=%.3f" % (a_obs, p))

print("\n=== D) per-excitation (shape-normalized, PCA8) ===")
for k, ex in enumerate(EXC):
    Xk = C[:, :, k].copy(); Xk = Xk / Xk.sum(1, keepdims=True)
    a = loo_auc(Xk, y, 8, "lr")
    print("  ex%3d  LOO AUC=%.3f" % (ex, a))

print("\n=== E) broad-dimming + emission peak-position shift ===")
tot = C.sum((1, 2))
print("  total intensity sham %.3g vs TBI %.3g  (MWU p=%.3g)" %
      (tot[y==0].mean(), tot[y==1].mean(), mannwhitneyu(tot[y==0], tot[y==1]).pvalue))
# spectral centroid (emission center of mass) per excitation
for k, ex in enumerate(EXC):
    sp = C[:, :, k]
    cen = (sp * em).sum(1) / sp.sum(1)
    d = cen[y==1].mean() - cen[y==0].mean()
    p = mannwhitneyu(cen[y==0], cen[y==1]).pvalue
    print("  ex%3d centroid sham %.2f TBI %.2f  shift %+.2f nm  (MWU p=%.3g)" %
          (ex, cen[y==0].mean(), cen[y==1].mean(), d, p))
