"""
exp3 analysis — organoid conditioned media, SHAM vs TBI (+ drug arms resveratrol / MG132),
two timepoints (3H, 6H) analysed separately (README objective 3).

Objective 1: discriminate TBI- from sham-conditioned media.
Objective 2: test drug effects by projecting drug arms onto the sham-vs-TBI axis
             (expect resveratrol -> rescue toward sham; MG132 -> worsen toward/beyond TBI).

Confound note: sham and TBI sit on SEPARATE plates, so condition is confounded with plate/batch.
Mitigations: (a) brightness-invariant SHAPE features (per-excitation L1 norm); (b) the drug arms
resveratrol/TBI/MG132 live on the SAME TBI plate, so if the cross-plate axis ALSO orders those
same-plate arms as predicted, that is biology not batch; (c) directional cosine test of group diffs.
"""
import os, sys, warnings
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[v] = "1"
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from repo_paths import data_path, results_dir, env_path
import numpy as np
import exp3_load as L
from scipy.stats import spearmanr
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def pr(*a): print(*a); sys.stdout.flush()

OUT = env_path("EXP3_OUT", "") or results_dir("04_drug_response_prediction"); os.makedirs(OUT, exist_ok=True)
EXC = L.EXC
COL = {"sham":"#2e86c1","tbi":"#cb4335","resv_tbi":"#27ae60","mg132_tbi":"#8e44ad",
       "resv_sham":"#a9cce3","mg132_sham":"#d2b4de","swcnt_noncm":"#7f8c8d"}

def shape_feats(cube, idx):
    """per-excitation L1-normalized shape spectra -> (n_wells, 512*5)."""
    X = cube[idx].astype(float)                       # (n,512,5)
    X = X / X.sum(1, keepdims=True)                   # L1 per (well,excitation)
    return X.reshape(len(idx), -1), X

def load_tp(tp):
    cs, em = L.load_plate("SHAM", tp); Gs = L.groups("SHAM")
    ct, _  = L.load_plate("TBI",  tp); Gt = L.groups("TBI")
    return em, cs, Gs, ct, Gt

def _loo_lda_auc(Z, y):
    """LOO-LDA AUC in a (small) feature space Z."""
    sc = np.zeros(len(y))
    lda = LDA()
    for tr, te in LeaveOneOut().split(Z):
        lda.fit(Z[tr], y[tr]); sc[te] = lda.decision_function(Z[te])
    return roc_auc_score(y, sc)

def loo_auc(Xa, Xb, ncomp=20):
    """Reduce to a fixed UNSUPERVISED PCA subspace (label-blind, so valid for permutation),
    then LOO-LDA. Same pipeline used for observed and permuted -> exact null."""
    X = np.vstack([Xa, Xb]); y = np.r_[np.zeros(len(Xa)), np.ones(len(Xb))]
    nc = min(ncomp, len(y) - 2)
    Z = make_pipeline_pca(X, nc)
    return _loo_lda_auc(Z, y), Z, y

def make_pipeline_pca(X, nc):
    sc = StandardScaler().fit(X); pca = PCA(nc, random_state=0).fit(sc.transform(X))
    return pca.transform(sc.transform(X))

def perm_null(Z, y, obs, n=2000, seed=0):
    rng = np.random.default_rng(seed); cnt = 0
    for _ in range(n):
        yp = rng.permutation(y)
        if _loo_lda_auc(Z, yp) >= obs: cnt += 1
    return (cnt + 1) / (n + 1)

def axis_project(Xsham, Xtbi, groups_feat):
    """fit sham-vs-TBI LDA axis on full data, oriented sham<low TBI>high; project groups.
    returns dict name->scores, plus the axis weight vector (in standardized space mapped back)."""
    X = np.vstack([Xsham, Xtbi]); y = np.r_[np.zeros(len(Xsham)), np.ones(len(Xtbi))]
    sc = StandardScaler().fit(X); lda = LDA().fit(sc.transform(X), y)
    def proj(M): return lda.decision_function(sc.transform(M))
    s0, s1 = proj(Xsham).mean(), proj(Xtbi).mean()
    flip = 1.0 if s1 > s0 else -1.0
    out = {"sham": flip*proj(Xsham), "tbi": flip*proj(Xtbi)}
    for n, M in groups_feat.items(): out[n] = flip*proj(M)
    return out

def cohend(a, b):
    na, nb = len(a), len(b); s = np.sqrt(((na-1)*a.var(ddof=1)+(nb-1)*b.var(ddof=1))/(na+nb-2))
    return (b.mean()-a.mean())/s

# ============================ run both timepoints ============================
results = {}
for tp in ("3H", "6H"):
    em, cs, Gs, ct, Gt = load_tp(tp)
    Xsham, Rsham = shape_feats(cs, Gs["sham"])
    Xtbi,  Rtbi  = shape_feats(ct, Gt["tbi"])
    # other arms
    feats = {}
    raw = {}
    for nm, (cube, idx) in {"resv_tbi":(ct,Gt["resv_tbi"]), "mg132_tbi":(ct,Gt["mg132_tbi"]),
                            "resv_sham":(cs,Gs["resv_sham"]), "mg132_sham":(cs,Gs["mg132_sham"]),
                            "swcnt_noncm":(cs,Gs["swcnt_noncm"])}.items():
        feats[nm], raw[nm] = shape_feats(cube, idx)

    pr("\n################  exp3  %s  ################" % tp)
    # ---- Objective 1: sham vs TBI ----
    auc, Z, y = loo_auc(Xsham, Xtbi)
    p = perm_null(Z, y, auc, n=2000)
    pr("[obj1] sham(n=%d) vs TBI(n=%d)  LOO-LDA AUC=%.3f  perm-p=%.4f" % (len(Xsham), len(Xtbi), auc, p))
    # mean shape diff (avg over excitation) and chirality of change
    mdiff = (Rtbi.mean(0) - Rsham.mean(0)).mean(1)    # (512,) avg over excitation
    gain = em[np.argsort(mdiff)[-5:]]; loss = em[np.argsort(mdiff)[:5]]
    cen_s = (em @ Rsham.mean(0).mean(1)) / Rsham.mean(0).mean(1).sum()
    cen_t = (em @ Rtbi.mean(0).mean(1)) / Rtbi.mean(0).mean(1).sum()
    pr("       centroid sham=%.1fnm TBI=%.1fnm  shift=%+.2fnm | gain@%s loss@%s"
       % (cen_s, cen_t, cen_t-cen_s, np.round(np.sort(gain)).astype(int), np.round(np.sort(loss)).astype(int)))

    # ---- Objective 2: drug projection onto TBI axis ----
    P = axis_project(Xsham, Xtbi, feats)
    pr("[obj2] projection on sham(0)->TBI(+) axis  (mean +/- sd):")
    for nm in ["sham","resv_sham","mg132_sham","tbi","resv_tbi","mg132_tbi","swcnt_noncm"]:
        v = P[nm]; pr("       %-12s % .3f +/- %.3f" % (nm, v.mean(), v.std()))
    # rescue / worsen effect sizes (same TBI plate -> batch-controlled)
    d_resv = cohend(P["tbi"], P["resv_tbi"])     # negative => toward sham (rescue)
    d_mg   = cohend(P["tbi"], P["mg132_tbi"])    # positive => beyond TBI (worse)
    pr("       resv_tbi vs tbi  Cohen d=%+.2f (neg=rescue)   mg132_tbi vs tbi d=%+.2f (pos=worse)"
       % (d_resv, d_mg))
    # directional cosine: does (resv_tbi - tbi) oppose (TBI - sham)?  (same-plate diff)
    vTBI = Xtbi.mean(0) - Xsham.mean(0)
    def cos(a, b): return float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30))
    cos_resv = cos(feats["resv_tbi"].mean(0)-Xtbi.mean(0), vTBI)
    cos_mg   = cos(feats["mg132_tbi"].mean(0)-Xtbi.mean(0), vTBI)
    pr("       cos((resv_tbi-tbi), (TBI-sham))=%+.2f (neg=rescues along axis) | cos(mg132)=%+.2f"
       % (cos_resv, cos_mg))

    results[tp] = dict(em=em, Rsham=Rsham, Rtbi=Rtbi, mdiff=mdiff, P=P, auc=auc, p=p,
                       cen_s=cen_s, cen_t=cen_t, d_resv=d_resv, d_mg=d_mg,
                       cos_resv=cos_resv, cos_mg=cos_mg, vTBI=vTBI)

# ---- maturation 3H vs 6H ----
def cos(a, b): return float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30))
cmat = cos(results["3H"]["mdiff"], results["6H"]["mdiff"])
pr("\n[maturation] cos(TBI-sham shape diff 3H, 6H) = %+.3f" % cmat)

# ============================ figures ============================
fig, axes = plt.subplots(2, 3, figsize=(17, 9))
for ri, tp in enumerate(("3H", "6H")):
    R = results[tp]; em = R["em"]
    # col0: mean shape spectra sham vs TBI (avg excitation)
    ax = axes[ri, 0]
    ax.plot(em, R["Rsham"].mean(0).mean(1), c=COL["sham"], lw=2, label="sham")
    ax.plot(em, R["Rtbi"].mean(0).mean(1),  c=COL["tbi"],  lw=2, label="TBI")
    ax.set_xlim(900, 1400); ax.set_xlabel("emission (nm)"); ax.set_ylabel("shape (L1-norm)")
    ax.set_title("[%s] mean shape  AUC=%.3f p=%.3f" % (tp, R["auc"], R["p"]), fontsize=10)
    ax.legend(fontsize=8)
    # col1: difference spectrum TBI-sham
    ax = axes[ri, 1]
    ax.axhline(0, c="k", lw=.6)
    ax.plot(em, R["mdiff"], c="#34495e", lw=1.5)
    ax.fill_between(em, R["mdiff"], 0, where=R["mdiff"]>0, color=COL["tbi"], alpha=.4)
    ax.fill_between(em, R["mdiff"], 0, where=R["mdiff"]<0, color=COL["sham"], alpha=.4)
    ax.set_xlim(900, 1400); ax.set_xlabel("emission (nm)"); ax.set_ylabel("TBI - sham shape")
    ax.set_title("[%s] diff: centroid %+.2fnm" % (tp, R["cen_t"]-R["cen_s"]), fontsize=10)
    # col2: drug projection strip
    ax = axes[ri, 2]; P = R["P"]
    order = ["sham","resv_sham","mg132_sham","tbi","resv_tbi","mg132_tbi","swcnt_noncm"]
    for xi, nm in enumerate(order):
        v = P[nm]
        ax.scatter(np.full(len(v), xi)+np.random.uniform(-.12,.12,len(v)), v,
                   c=COL[nm], s=26, edgecolor="k", lw=.3, zorder=3)
        ax.plot([xi-.25, xi+.25], [v.mean(), v.mean()], c="k", lw=2, zorder=4)
    ax.axhline(P["sham"].mean(), c=COL["sham"], ls="--", lw=1)
    ax.axhline(P["tbi"].mean(),  c=COL["tbi"],  ls="--", lw=1)
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=40, ha="right", fontsize=8)
    ax.set_ylabel("sham(0) -> TBI(+) axis score")
    ax.set_title("[%s] drug arms  resvD=%+.2f mgD=%+.2f" % (tp, R["d_resv"], R["d_mg"]), fontsize=10)
fig.suptitle("exp3 organoid conditioned media — SHAM vs TBI + drug arms (3H vs 6H separate)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
p1 = os.path.join(OUT, "exp3_sham_tbi_drugs.png"); fig.savefig(p1, dpi=150)
pr("\nsaved", p1)
pr("[maturation] cos 3H/6H diff spectra = %+.3f" % cmat)
