"""Two polished figures (corrected labels: col10=1/4 TBI, col12=3/4 TBI):
  A) sham-vs-TBI LDA decision-score dose-response across negctrl/sham/mixes/TBI
  B) PCA of shape spectra, all wells, showing mixes nesting between sham and TBI
Run with MIX_INVERTED=1.
"""
import os
os.environ.setdefault("MIX_INVERTED", "1")
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
TAG = os.environ.get("EXP2_TAG", "2h")
import numpy as np, sys
np.random.seed(0)  # legacy np.random jitter in figures
from exp2_load import load_all, labels
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.metrics import roc_auc_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
frac   = {"negctrl":np.nan,"sham":0.0,"mix25":0.25,"mix50":0.5,"mix75":0.75,"TBI":1.0}
labtxt = {"negctrl":"negctrl","sham":"sham (0%)","mix25":"¼ TBI","mix50":"½ TBI","mix75":"¾ TBI","TBI":"TBI (100%)"}
colmap = {"negctrl":"#7f8c8d","sham":"#2e86c1","mix25":"#aed6a1","mix50":"#5dade2",
          "mix75":"#239b56","TBI":"#cb4335"}
# dose green gradient for mixes
colmap.update({"sham":"#2e86c1","mix25":"#a9dfbf","mix50":"#52be80","mix75":"#1e8449","TBI":"#cb4335"})

Xall = {g: shp(cube[wells(g)]) for g in groups}

# fit PCA on sham+TBI shape, then LDA decision axis
tr = np.concatenate([wells("sham"), wells("TBI")])
ytr = np.r_[np.zeros(len(wells("sham"))), np.ones(len(wells("TBI")))]
pca = make_pipeline(StandardScaler(), PCA(20, random_state=0)).fit(shp(cube[tr]))
Ztr = pca.transform(shp(cube[tr]))
lda = LDA().fit(Ztr, ytr)
def dscore(g): return lda.decision_function(pca.transform(Xall[g]))

# ================= FIGURE A: decision-score dose-response =================
figA, ax = plt.subplots(figsize=(7.6, 5.4))
xs = [frac[g] for g in ["sham","mix25","mix50","mix75","TBI"]]
means = [];
for g in ["sham","mix25","mix50","mix75","TBI"]:
    s = dscore(g); means.append(s.mean())
    ax.scatter(np.full(len(s), frac[g]) + np.random.uniform(-.012,.012,len(s)),
               s, c=colmap[g], s=46, edgecolor="k", lw=.4, zorder=3)
ax.plot(xs, means, "-", c="#34495e", lw=2.2, zorder=2)
ax.scatter(xs, means, c="k", s=70, marker="D", zorder=4, label="group mean")
# negctrl reference band
nc = dscore("negctrl")
ax.axhspan(nc.mean()-nc.std(), nc.mean()+nc.std(), color="#7f8c8d", alpha=.15, zorder=0)
ax.axhline(nc.mean(), c="#7f8c8d", ls=":", lw=1.6, label="negctrl (SWCNT in media)")
ax.axhline(0, c="grey", lw=.6)
ax.set_xlabel("TBI-media fraction", fontsize=12)
ax.set_ylabel("sham-vs-TBI LDA decision score\n(← sham   |   TBI →)", fontsize=12)
ax.set_xticks(xs); ax.set_xticklabels(["0\nsham","¼","½","¾","1\nTBI"])
ax.set_title("exp2 / %s — dilution dose-response\n(classifier trained on sham & TBI only)" % TAG, fontsize=12)
ax.legend(loc="upper left", fontsize=9)
figA.tight_layout()
OUT = os.environ.get("EXP2_OUT", "/mnt/c/Users/riccardo-s/Documents/CNT/targetALS/exp5/tbi_sham_tbi_outputs")
os.makedirs(OUT, exist_ok=True)
pA = os.path.join(OUT, "exp2_%s_decision_dose_response.png" % TAG); figA.savefig(pA, dpi=160)
pr("saved", pA)

# ================= FIGURE B: PCA (2 panels) =================
from matplotlib.patches import Ellipse
def ellipse(ax, X, color):
    if len(X) < 3: return
    mu = X.mean(0); cov = np.cov(X.T); val, vec = np.linalg.eigh(cov)
    ang = np.degrees(np.arctan2(*vec[:, -1][::-1]))
    w, h = 2 * np.sqrt(np.maximum(val, 0))
    ax.add_patch(Ellipse(mu, w, h, angle=ang, facecolor=color, alpha=.12, edgecolor=color, lw=1.5, zorder=1))
order_plot = ["negctrl","sham","mix25","mix50","mix75","TBI"]
traj = ["sham","mix25","mix50","mix75","TBI"]

# --- panel 1: full-spectrum PCA (variance PC0 vs PC1) — honest, multivariate ---
Zfull = {g: pca.transform(Xall[g]) for g in groups}
# --- panel 2: driver-band PCA (emission 1080-1200 nm) — PC0 = sham->TBI dose axis ---
band = (em >= 1080) & (em <= 1200)
def shp_band(C):
    S = C.copy()
    for k in range(C.shape[2]): S[:, :, k] /= S[:, :, k].sum(1, keepdims=True)
    return S[:, band, :].reshape(len(C), -1)
scb = StandardScaler().fit(shp_band(cube[tr])); pcb = PCA(5, random_state=0).fit(scb.transform(shp_band(cube[tr])))
Zb = {g: pcb.transform(scb.transform(shp_band(cube[wells(g)]))) for g in groups}
# orient PC0 so TBI is positive
if Zb["TBI"][:,0].mean() < Zb["sham"][:,0].mean():
    for g in groups: Zb[g][:,0] *= -1
shb, tbb = Zb["sham"], Zb["TBI"]
auc0 = roc_auc_score(np.r_[np.zeros(len(shb)),np.ones(len(tbb))], np.r_[shb[:,0],tbb[:,0]])
evb = pcb.explained_variance_ratio_

figB, axes = plt.subplots(1, 2, figsize=(14, 6))
for ax, (Zd, p1, p2, ttl, xl, yl) in zip(axes, [
        (Zfull, 0, 1, "Full shape-spectrum PCA\n(no 2 PCs separate — signal is multivariate; full 20-PC classifier AUC=1.0)",
         "PC1 (%.0f%% var)" % (pca.named_steps['pca'].explained_variance_ratio_[0]*100),
         "PC2 (%.0f%% var)" % (pca.named_steps['pca'].explained_variance_ratio_[1]*100)),
        (Zb, 0, 1, "PCA of mid-diameter driver band (1080–1200 nm)\nPC1 = sham→TBI dose axis (AUC %.2f)" % auc0,
         "PC1  (driver band, %.0f%% var)" % (evb[0]*100), "PC2  (%.0f%% var)" % (evb[1]*100))]):
    for g in order_plot:
        Z = Zd[g]
        ax.scatter(Z[:,p1], Z[:,p2], c=colmap[g], s=64, edgecolor="k", lw=.5, alpha=.9, zorder=3, label=labtxt[g])
    ellipse(ax, Zd["sham"][:, [p1,p2]], "#2e86c1"); ellipse(ax, Zd["TBI"][:, [p1,p2]], "#cb4335")
    cents = {g: Zd[g].mean(0) for g in order_plot}
    tx = [cents[g][p1] for g in traj]; ty = [cents[g][p2] for g in traj]
    ax.plot(tx, ty, "-", c="k", lw=1.8, zorder=4)
    for g in traj:
        ax.scatter(cents[g][p1], cents[g][p2], facecolor="none", edgecolor="k", s=230, lw=2, zorder=5)
    ax.set_xlabel(xl, fontsize=11); ax.set_ylabel(yl, fontsize=11); ax.set_title(ttl, fontsize=11)
axes[0].legend(fontsize=8.5, loc="best", framealpha=.9)
figB.suptitle("exp2 / %s — shape-spectrum PCA  (dilution mixes ¼→¾ in green, dose path = black)" % TAG, fontsize=13)
figB.tight_layout(rect=[0,0,1,0.96])
pB = os.path.join(OUT, "exp2_%s_pca_dose.png" % TAG); figB.savefig(pB, dpi=160)
pr("saved", pB, "| driver-band PC1 AUC %.2f, centroids(sham..TBI): %s" %
   (auc0, " ".join("%+.1f"%Zb[g][:,0].mean() for g in traj)))
