"""Compare the dilution dose-response across classifiers (incl. elastic-net).
All trained on sham(32) vs TBI(32) shape-normalized spectra; mixes/negctrl projected.
Each classifier's score is rescaled so sham mean=0, TBI mean=1 -> common dose axis.
Run with MIX_INVERTED=1 (col10=1/4 TBI, col12=3/4 TBI).
"""
import os
os.environ.setdefault("MIX_INVERTED", "1")
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
TAG = os.environ.get("EXP2_TAG", "2h")
import numpy as np, sys
from exp2_load import load_all, labels
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.model_selection import LeaveOneOut
from sklearn.base import clone
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr, kendalltau
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def pr(*a): print(*a); sys.stdout.flush()

cube, em = load_all(); rows, cols, grp = labels()
def shp(C):
    S = C.copy()
    for k in range(C.shape[2]): S[:, :, k] /= S[:, :, k].sum(1, keepdims=True)
    return S.reshape(len(C), -1)
def wells(g):
    idx = np.where(grp == g)[0]
    return idx[idx != 0] if g == "negctrl" else idx

groups = ["negctrl","sham","mix25","mix50","mix75","TBI"]
order  = ["sham","mix25","mix50","mix75","TBI"]
frac   = {"negctrl":np.nan,"sham":0.0,"mix25":0.25,"mix50":0.5,"mix75":0.75,"TBI":1.0}
X = {g: shp(cube[wells(g)]) for g in groups}
tr = np.concatenate([wells("sham"), wells("TBI")])
Xtr = shp(cube[tr]); ytr = np.r_[np.zeros(len(wells("sham"))), np.ones(len(wells("TBI")))]

def raw_score(model, has_df=True):
    """fit on train, return dict of raw decision scores per group + LOO AUC."""
    model.fit(Xtr, ytr)
    sc = {}
    for g in groups:
        sc[g] = model.decision_function(X[g]) if has_df else model.predict_proba(X[g])[:,1]
    # LOO AUC on train
    pr_ = np.zeros(len(ytr))
    for trn, te in LeaveOneOut().split(Xtr):
        m = clone(model).fit(Xtr[trn], ytr[trn])
        pr_[te] = (m.decision_function(Xtr[te]) if has_df else m.predict_proba(Xtr[te])[:,1])
    return sc, roc_auc_score(ytr, pr_)

# classifier pipelines (all on shape features)
PCA20 = lambda est: make_pipeline(StandardScaler(), PCA(20, random_state=0), est)
clfs = {
    "LDA (PCA20)":            (PCA20(LDA()), True),
    "LogReg L2 (PCA20)":      (PCA20(LogisticRegression(max_iter=5000, C=1.0)), True),
    "Linear SVM (PCA20)":     (PCA20(SVC(kernel="linear", C=0.5, random_state=0)), True),
    "Elastic-net (PCA20)":    (PCA20(LogisticRegression(penalty="elasticnet", solver="saga", random_state=0,
                                     l1_ratio=0.5, C=1.0, max_iter=20000)), True),
    "Elastic-net (full 2560)":(make_pipeline(StandardScaler(),
                                  LogisticRegression(penalty="elasticnet", solver="saga", random_state=0,
                                     l1_ratio=0.5, C=0.5, max_iter=20000)), True),
}

results = {}
pr("%-24s %6s | %-32s | %8s %8s" % ("classifier","LOO_AUC","group means (sham..TBI, rescaled 0..1)","Sp_all","Sp_mix"))
pr("-"*100)
for name, (model, has_df) in clfs.items():
    sc, auc = raw_score(model, has_df)
    s0, s1 = sc["sham"].mean(), sc["TBI"].mean()
    resc = {g: (sc[g] - s0) / (s1 - s0) for g in groups}      # sham->0, TBI->1
    means = {g: resc[g].mean() for g in groups}
    fr_all = np.concatenate([np.full(len(wells(g)), frac[g]) for g in order])
    v_all  = np.concatenate([resc[g] for g in order])
    sp_all = spearmanr(fr_all, v_all)[0]
    frm = np.concatenate([np.full(len(wells(g)), frac[g]) for g in ["mix25","mix50","mix75"]])
    vm  = np.concatenate([resc[g] for g in ["mix25","mix50","mix75"]])
    sp_mix = spearmanr(frm, vm)[0]
    results[name] = dict(auc=auc, resc=resc, means=means, sp_all=sp_all, sp_mix=sp_mix,
                         kt_all=kendalltau(fr_all, v_all)[0])
    gm = " ".join("%+.2f"%means[g] for g in order)
    pr("%-24s %6.3f | %-32s | %+7.2f %+7.2f" % (name, auc, gm, sp_all, sp_mix))

# cross-classifier agreement of per-well scores
pr("\nPearson correlation of rescaled per-well scores between classifiers (all 96 wells):")
names = list(clfs)
allv = {n: np.concatenate([results[n]["resc"][g] for g in groups]) for n in names}
for i in range(len(names)):
    row = " ".join("%.2f"%np.corrcoef(allv[names[i]], allv[names[j]])[0,1] for j in range(len(names)))
    pr("  %-24s %s" % (names[i], row))

# ---- figure: overlay dose-response per classifier ----
fig, ax = plt.subplots(figsize=(9, 6))
xs = [frac[g] for g in order]
palette = {"LDA (PCA20)":"#cb4335","LogReg L2 (PCA20)":"#2e86c1","Linear SVM (PCA20)":"#16a085",
           "Elastic-net (PCA20)":"#8e44ad","Elastic-net (full 2560)":"#e67e22"}
for name in names:
    r = results[name]; ms = [r["means"][g] for g in order]
    sd = [r["resc"][g].std()/np.sqrt(len(r["resc"][g])) for g in order]
    ax.errorbar(xs, ms, yerr=sd, marker="o", lw=2, capsize=3, color=palette[name],
                label="%s  (AUC %.2f, ρmix %+.2f)" % (name, r["auc"], r["sp_mix"]))
ax.plot([0,1],[0,1], ls=":", c="grey", lw=1.4, label="ideal linear dose")
ax.axhline(0, c="grey", lw=.5); ax.axhline(1, c="grey", lw=.5)
ax.set_xticks(xs); ax.set_xticklabels(["0\nsham","¼","½","¾","1\nTBI"])
ax.set_xlabel("TBI-media fraction", fontsize=12)
ax.set_ylabel("rescaled classifier score  (sham=0, TBI=1)", fontsize=12)
ax.set_title("exp2 / %s — dilution dose-response across classifiers\n(trained sham vs TBI; mixes projected; mean ± SEM)" % TAG, fontsize=12)
ax.legend(fontsize=8.5, loc="upper left")
fig.tight_layout()
OUT = os.environ.get("EXP2_OUT", "/mnt/c/Users/riccardo-s/Documents/CNT/targetALS/exp5/tbi_sham_tbi_outputs")
os.makedirs(OUT, exist_ok=True)
p = os.path.join(OUT, "exp2_%s_classifier_dose_compare.png" % TAG); fig.savefig(p, dpi=160)
pr("\nsaved", p)
