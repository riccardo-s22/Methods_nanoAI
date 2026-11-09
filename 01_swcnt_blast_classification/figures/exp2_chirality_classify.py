"""
exp2 (blast-conditioned-media SWCNT experiment) -- CHIRALITY peak-feature classification,
mirroring exactly what was done for exp3-3H, so the two experiments are directly comparable.

WHAT IT DOES
  For the 2 h and the 24 h exp2 reads, using ONLY the pre-computed split-Gaussian peak caches:
    1. sham vs blast(TBI) leave-one-out LDA AUC on the 66 pairwise chirality log10 intensity
       ratios built from the fitted Gaussian MAXIMA (MX), reduced to 15 PCs.
       Pipeline = StandardScaler -> PCA -> LDA, refitted INSIDE every LOO fold (no leakage).
    2. same for the 12 single-chirality AUC-fractions (PCA 8; 12 features cannot give 15 PCs)
       and for the combined 78-feature set (PCA 15).
    3. label-permutation test of the 66-ratio LOO AUC (B = 5000, fixed seed).
    4. dilution dose-response: the 66-ratio model refitted on ALL sham+blast wells, then the
       25/50/75 % blast dilution wells and the SWCNT-in-media negative control scored WITHOUT
       retraining; scores rescaled so sham mean = 0 and blast mean = 1.
    5. per-chirality and per-ratio sham->blast Cohen's d, with per-feature AUC, t-test p and
       BH-FDR q (same table layout as TBI/exp3/exp3_3H_peakfit_regularized_results.txt).

  Reference numbers it is being compared against (not recomputed here):
    exp2 full "shape" vectors (2560 feats) LOO AUC = 1.00 at 2 h and 24 h
    exp3-3H: 12 fractions 0.613 | 66 ratios 0.938 | 78 combined 0.920

CACHES READ (nothing is re-fitted for steps 1-3 and 5)
  /mnt/c/Users/riccardo-s/Documents/CNT/targetALS/exp5/exp2_2h_peakfeats.npz
  /mnt/c/Users/riccardo-s/Documents/CNT/targetALS/exp5/exp2_24h_peakfeats.npz
  keys: MXs/AUs = sham (32 wells x 12 chiralities), MXt/AUt = blast/TBI (32 x 12).
  Chirality column order is list(CHIRALITY_POSITIONS.keys()) from
  targetALS/exp5/physical_descriptors/chirality_interpolated_extraction.py, i.e.
  (8,3) (6,5) (7,5) (10,2) (9,4) (8,4) (7,6) (8,6) (8,7) (9,5) (10,3) (10,5).

  *** DEVIATION, step 4 only ***
  The exp2 caches contain ONLY sham and blast. They hold NO dilution and NO negative-control
  wells, so the dose-response cannot be done from cache alone. This script therefore extracts
  the split-Gaussian maxima for the 24 dilution wells + 8 negctrl wells per timepoint with the
  canonical extractor and stores them in exp2_<tp>_mix_peakfeats.npz next to this file (built
  once, reused afterwards). Validation: re-extracting the 64 CACHED sham/blast wells with the
  same code in this environment reproduces the cached MX to a median relative error of ~1e-11
  (95th pct ~1e-3), so the dilution MX values are on the same footing as the cached ones.
  The sham/blast training features themselves are always the CACHED ones, never re-fitted.
  (AU does not reproduce as cleanly -- a handful of wells land in a different local optimum --
  which is why only MX-based ratios, never AU fractions, are used for the dilution analysis.)

OUTPUTS (written next to this script)
  exp2_chirality_results.txt
  exp2_chirality.png / exp2_chirality.svg

DILUTION PLATE LABELLING: exp2_load.py defaults to col10 = 75 % blast / col12 = 25 % blast but
  notes the two may be inverted; every exp2 analysis in the paper (exp2_dose_pca_figs.py,
  exp2_classifier_compare.py, exp2_unsupervised*.py, paper/svg/README.txt, the concept figure)
  runs with MIX_INVERTED=1, i.e. col10 = 25 %, col11 = 50 %, col12 = 75 % blast. This script
  sets MIX_INVERTED=1 too, so the dose-response is directly comparable with the published
  shape-vector dose-response. The opposite labelling is reported as a sensitivity.

SEED: 0 everywhere (np.random.default_rng(0) for permutations, random_state=0 for PCA).
RUN:  /home/riccardo-s/miniforge3/bin/python exp2_chirality_classify.py
"""
import os, sys, warnings, importlib.util
from itertools import combinations
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"
warnings.filterwarnings("ignore")
import numpy as np
from scipy.stats import ttest_ind, spearmanr, sem
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score
from joblib import Parallel, delayed
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED = 0
B_PERM = 5000
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = "/mnt/c/Users/riccardo-s/Documents/CNT/targetALS/exp5"
EXP2DIR = "/mnt/c/Users/riccardo-s/Documents/CNT/TBI/exp2"
PD = "/mnt/c/Users/riccardo-s/Documents/CNT/targetALS/exp5/physical_descriptors"

# canonical chirality column order (= list(cie.CHIRALITY_POSITIONS.keys()))
ORDER = ["ch8_3", "ch6_5", "ch7_5", "ch10_2", "ch9_4", "ch8_4",
         "ch7_6", "ch8_6", "ch8_7", "ch9_5", "ch10_3", "ch10_5"]
EM_NM = [973.98, 987.82, 1047.81, 1080.60, 1131.96, 1130.34,
         1138.19, 1200.03, 1288.27, 1262.98, 1267.70, 1282.97]
SH = {k: k.replace("ch", "(").replace("_", ",") + ")" for k in ORDER}
CH = [SH[k] for k in ORDER]
PAIRS = list(combinations(range(12), 2))                      # 66 pairs
PNAME = ["%s/%s" % (CH[a], CH[b]) for a, b in PAIRS]
FEAT78 = ["frac%s" % c for c in CH] + ["r:%s" % n for n in PNAME]
flo = 1e-6

LOG = []
def pr(*a):
    s = " ".join(str(x) for x in a)
    print(s); sys.stdout.flush(); LOG.append(s)

# ---------------------------------------------------------------- features
def ratios(MX):
    """66 pairwise log10 intensity ratios -- exactly TBI/exp3/exp3_drug_axis.py lines 23-26."""
    return np.column_stack([np.log10((MX[:, a] + flo) / (MX[:, b] + flo)) for a, b in PAIRS])

def fracs(AU):
    """12 single-chirality AUC fractions -- exactly TBI/exp3/exp3_peakfit.py frac()."""
    A = np.nan_to_num(np.asarray(AU, float), nan=0.0); A[A < 0] = 0
    return A / np.clip(A.sum(1, keepdims=True), 1e-30, None)

# ---------------------------------------------------------------- CV helpers
def loo_scores(X, y, npc):
    """LOO decision_function scores; scaler+PCA+LDA refitted inside every fold."""
    sc = np.zeros(len(y))
    for tr, te in LeaveOneOut().split(X):
        pipe = make_pipeline(StandardScaler(), PCA(n_components=npc, random_state=SEED), LDA())
        pipe.fit(X[tr], y[tr])
        sc[te] = pipe.decision_function(X[te])
    return sc

def loo_auc(X, y, npc):
    return roc_auc_score(y, loo_scores(X, y, npc))

def cohend(a, b):
    s = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return (b.mean() - a.mean()) / (s + 1e-12)

def bh(p):
    p = np.asarray(p, float); o = np.argsort(p)
    ranked = p[o] * len(p) / (np.arange(len(p)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(len(p)); q[o] = np.clip(ranked, 0, 1); return q

# ---------------------------------------------------------------- dilution MX (step 4 only)
def mix_MX(tp):
    """MX for the 25/50/75 % dilution wells + negctrl. Cached in this folder after first build."""
    cpath = os.path.join(HERE, "exp2_%s_mix_peakfeats.npz" % tp)
    if os.path.exists(cpath):
        z = np.load(cpath)
        return {k: z[k] for k in z.files}
    if not hasattr(np, "trapz"):
        np.trapz = np.trapezoid                               # numpy>=2 shim for the extractor
    spec = importlib.util.spec_from_file_location(
        "cie", os.path.join(PD, "chirality_interpolated_extraction.py"))
    cie = importlib.util.module_from_spec(spec); spec.loader.exec_module(cie)
    CHIR = cie.CHIRALITY_POSITIONS
    EXC = np.array([570, 640, 670, 750, 780], float); exmin, exmax = EXC.min(), EXC.max()
    os.environ.setdefault("MIX_INVERTED", "1")                # paper convention, see below
    sys.path.insert(0, EXP2DIR)
    import exp2_load
    exp2_load.DIR = os.path.join(EXP2DIR, tp)
    cube, em = exp2_load.load_all()
    rows, cols, grp = exp2_load.labels()

    def well(mat):
        sp, emf, exf = cie.build_interpolated_surface(mat, em, EXC, factor=cie.INTERP_FACTOR)
        mx = np.zeros(12)
        for k, (name, c) in enumerate(CHIR.items()):
            tex = min(max(c["excitation_nm"], exmin), exmax)
            be, bx, m = cie.find_local_max_in_radius(
                sp, emf, exf, c["emission_nm"], tex, cie.SEARCH_RADIUS_NM)
            ep, ip = cie.extract_consistent_gradient_points(sp, emf, be, bx)
            mx[k] = max(cie.fit_gaussian_and_extract(ep, ip, be, m)["gauss_max"], 0)
        return mx

    out = {}
    for g in ["mix25", "mix50", "mix75", "negctrl"]:
        idx = [i for i in np.where(grp == g)[0] if i != 0]    # drop A1 ramp placeholder
        out[g] = np.array([well(cube[w].astype(float)) for w in idx])
    np.savez(cpath, **out)
    return out

# ---------------------------------------------------------------- per timepoint
def run(tp):
    pr("\n" + "#" * 78)
    pr("####  exp2 %s -- chirality peak-feature classification" % tp.upper())
    pr("#" * 78)
    z = np.load(os.path.join(CACHE, "exp2_%s_peakfeats.npz" % tp))
    pr("  cache: exp2_%s_peakfeats.npz  keys/shapes: %s" % (
        tp, ", ".join("%s%s" % (k, z[k].shape) for k in z.files)))
    MXs, AUs, MXt, AUt = z["MXs"], z["AUs"], z["MXt"], z["AUt"]
    ns, nt = len(MXs), len(MXt)
    pr("  sham n=%d   blast(TBI) n=%d   (A1 acquisition placeholder already absent)" % (ns, nt))

    y = np.r_[np.zeros(ns), np.ones(nt)]
    Rs, Rt = ratios(MXs), ratios(MXt)
    Fs, Ft = fracs(AUs), fracs(AUt)
    XR = np.vstack([Rs, Rt])                                  # 66
    XF = np.vstack([Fs, Ft])                                  # 12
    XB = np.hstack([np.vstack([Fs, Ft]), np.vstack([Rs, Rt])])  # 78

    # ---- 1/2 LOO AUCs
    pr("\n=== [1-2] leave-one-out LDA AUC (StandardScaler -> PCA -> LDA, refit per fold) ===")
    aR = loo_auc(XR, y, 15)
    aF = loo_auc(XF, y, 8)
    aB = loo_auc(XB, y, 15)
    pr("  [66 log-ratios  -> PCA15]  LOO-LDA AUC = %.3f" % aR)
    pr("  [12 fractions   -> PCA8 ]  LOO-LDA AUC = %.3f" % aF)
    pr("       (12 features cannot give 15 PCs; 8 chosen per spec. sensitivity:")
    for n in (4, 6, 10, 11, 12):
        pr("        PCA%-2d AUC = %.3f%s" % (n, loo_auc(XF, y, n),
            "   <- DEGENERATE, see note" if n == 12 else ""))
    evr = PCA(12).fit(StandardScaler().fit_transform(XF)).explained_variance_ratio_
    pr("        note: the 12 fractions are compositional (rows sum to 1) so the centred matrix")
    pr("              has rank 11; PC12 carries explained-variance ratio %.1e, i.e. machine" % evr[-1])
    pr("              precision noise. Any AUC gain at PCA12 is a numerical artefact of that")
    pr("              null direction and must not be quoted. PCA<=11 is the honest range.")
    pr("       )")
    pr("  [78 combined    -> PCA15]  LOO-LDA AUC = %.3f" % aB)
    pr("  reference: exp2 full shape vector (2560 feats) LOO AUC = 1.00 at this timepoint")
    pr("  reference: exp3-3H  fractions 0.613 | ratios 0.938 | combined 0.920")

    # ---- 3 permutation test on the 66 ratios
    pr("\n=== [3] label-permutation test, 66 log-ratios -> PCA15, B=%d, seed=%d ===" % (B_PERM, SEED))
    npath = os.path.join(HERE, "exp2_%s_permnull.npz" % tp)
    if os.path.exists(npath):
        null = np.load(npath)["null"]
        pr("  (null distribution reloaded from exp2_%s_permnull.npz)" % tp)
    else:
        rng = np.random.default_rng(SEED)
        YP = [rng.permutation(y) for _ in range(B_PERM)]       # drawn once, reproducibly
        def one(yp):
            return roc_auc_score(yp, loo_scores(XR, yp, 15))
        null = np.array(Parallel(n_jobs=min(8, os.cpu_count() or 1), batch_size=16)(
            delayed(one)(yp) for yp in YP))
        np.savez(npath, null=null)
    k = int((null >= aR).sum())
    pmc = (k + 1) / (B_PERM + 1)
    zsc = (aR - null.mean()) / (null.std(ddof=1) + 1e-12)
    pr("  observed AUC          = %.4f" % aR)
    pr("  null mean +/- sd      = %.4f +/- %.4f" % (null.mean(), null.std(ddof=1)))
    pr("  null median           = %.4f" % np.median(null))
    pr("  null 95th percentile  = %.4f    null max = %.4f" % (np.percentile(null, 95), null.max()))
    pr("  exceedances k (null >= obs) = %d / %d" % (k, B_PERM))
    pr("  Monte-Carlo p = (k+1)/(B+1) = %.2e" % pmc)
    pr("  z vs null     = %+.2f" % zsc)

    # ---- 4 dilution dose-response on the 66-ratio axis
    pr("\n=== [4] dilution dose-response, 66-ratio model fitted on ALL sham+blast, dilutions "
       "scored without retraining ===")
    pr("  NOTE: dilution/negctrl wells are NOT in the peak-feature cache (cache = sham+blast")
    pr("        only). Their MX was extracted once with the canonical extractor; the cached")
    pr("        sham/blast MX is used unchanged for training. See the module docstring.")
    MIX = mix_MX(tp)
    pipe = make_pipeline(StandardScaler(), PCA(n_components=15, random_state=SEED), LDA())
    pipe.fit(XR, y)
    s_sham = pipe.decision_function(Rs); s_blast = pipe.decision_function(Rt)
    m0, m1 = s_sham.mean(), s_blast.mean()
    resc = lambda s: (s - m0) / (m1 - m0)
    S = {"sham": resc(s_sham), "blast": resc(s_blast)}
    for g in ["mix25", "mix50", "mix75", "negctrl"]:
        if len(MIX[g]) == 0:
            pr("  !! group %s has 0 wells -- skipped" % g); continue
        S[g] = resc(pipe.decision_function(ratios(MIX[g])))
    pr("  (apparent in-sample separation of the training fit: sham %+.2f, blast %+.2f rescaled)"
       % (S["sham"].mean(), S["blast"].mean()))
    pr("\n  group        n   rescaled score (sham=0, blast=1)")
    for g in ["negctrl", "sham", "mix25", "mix50", "mix75", "blast"]:
        pr("   %-9s %3d   %+.3f +/- %.3f (SEM)" % (g, len(S[g]), S[g].mean(), sem(S[g])))
    BF = {"mix25": 0.25, "mix50": 0.50, "mix75": 0.75}
    xm = np.concatenate([[BF[g]] * len(S[g]) for g in ["mix25", "mix50", "mix75"]])
    ym = np.concatenate([S[g] for g in ["mix25", "mix50", "mix75"]])
    rho_mix, p_mix = spearmanr(xm, ym)
    xa = np.concatenate([[0.0] * len(S["sham"]), xm, [1.0] * len(S["blast"])])
    ya = np.concatenate([S["sham"], ym, S["blast"]])
    rho_all, p_all = spearmanr(xa, ya)
    pr("\n  Spearman rho_mix  (dilution wells only, n=%d) = %+.3f  p=%.2e" % (len(ym), rho_mix, p_mix))
    pr("  Spearman rho_all  (sham+dilutions+blast, n=%d) = %+.3f  p=%.2e" % (len(ya), rho_all, p_all))
    pr("  reference: shape-vector version at 24 h gave group means ~0.19/0.48/0.78 and")
    pr("             rho_mix 0.86-0.94 across five models")
    pr("  labelling: MIX_INVERTED=1 (col10=25%, col11=50%, col12=75% blast) -- the convention")
    pr("             used by every exp2 analysis in the paper. Sensitivity under the opposite")
    xm_inv = np.concatenate([[BF[g]] * len(S[g]) for g in ["mix75", "mix50", "mix25"]])
    pr("             (exp2_load.py default) labelling: rho_mix = %+.3f" % spearmanr(xm_inv, ym)[0])

    # ---- 5 effect sizes
    pr("\n=== [5] sham->blast effect sizes (+ = higher in blast) ===")
    pr("\n  per-chirality AUC-fraction:")
    pr("   %-8s %-9s %8s %8s %8s %8s %6s" % ("(n,m)", "em(nm)", "sham", "blast", "d", "AUC", "q"))
    pf, df_, af = [], [], []
    for k_ in range(12):
        a, b = Fs[:, k_], Ft[:, k_]
        df_.append(cohend(a, b)); pf.append(ttest_ind(a, b).pvalue)
        af.append(roc_auc_score(y, np.r_[a, b])); af[-1] = max(af[-1], 1 - af[-1])
    qf = bh(pf)
    for k_ in np.argsort(-np.abs(df_)):
        pr("   %-8s %-9.1f %8.4f %8.4f %+8.2f %8.2f %6.2f%s" % (
            CH[k_], EM_NM[k_], Fs[:, k_].mean(), Ft[:, k_].mean(),
            df_[k_], af[k_], qf[k_], " *" if qf[k_] < 0.05 else ""))
    pr("   # fractions with q<0.05: %d / 12" % int((qf < 0.05).sum()))

    dr, prr, ar_ = [], [], []
    for j in range(66):
        a, b = Rs[:, j], Rt[:, j]
        dr.append(cohend(a, b)); prr.append(ttest_ind(a, b).pvalue)
        v = roc_auc_score(y, np.r_[a, b]); ar_.append(max(v, 1 - v))
    dr = np.array(dr); prr = np.array(prr); ar_ = np.array(ar_); qr = bh(prr)
    pr("\n  top discriminative intensity ratios (sham vs blast), BH-FDR over 66:")
    pr("   # ratios with q<0.05: %d / 66" % int((qr < 0.05).sum()))
    top = np.argsort(-np.abs(dr))[:12]
    for j in top:
        pr("   %-16s d=%+.2f AUC=%.2f p=%.1e q=%.2e%s" % (
            PNAME[j], dr[j], ar_[j], prr[j], qr[j], " *" if qr[j] < 0.05 else ""))
    pr("\n  ratio |d| summary: median %.2f  mean %.2f  max %.2f  (# |d|>0.8: %d/66)" % (
        np.median(np.abs(dr)), np.abs(dr).mean(), np.abs(dr).max(), int((np.abs(dr) > 0.8).sum())))

    return dict(tp=tp, aucR=aR, aucF=aF, aucB=aB, null=null, k=k, p=pmc, z=zsc,
                S=S, dr=dr, ar=ar_, qr=qr, top=top, rho_mix=rho_mix, rho_all=rho_all,
                df=np.array(df_), qf=qf)

# ---------------------------------------------------------------- figure
def figure(res):
    plt.rcParams.update({"svg.fonttype": "none", "font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 8,
                         "axes.linewidth": 0.8, "pdf.fonttype": 42})
    fig, AX = plt.subplots(2, 3, figsize=(11.0, 6.4))
    CB = "#2e86c1"; CR = "#cb4335"
    for i, r in enumerate(res):
        tp = r["tp"]
        # (a) top ratios by |Cohen d|
        ax = AX[i, 0]
        top = np.argsort(-np.abs(r["dr"]))[:12][::-1]
        vals = r["dr"][top]
        ax.barh(np.arange(12), vals, color=[CR if v > 0 else CB for v in vals],
                height=0.72, edgecolor="none")
        ax.set_yticks(np.arange(12)); ax.set_yticklabels([PNAME[j] for j in top], fontsize=6.5)
        ax.axvline(0, c="k", lw=0.6)
        ax.set_xlabel("Cohen's $d$ (sham $\\rightarrow$ blast)")
        ax.set_title("%s  (a) top chirality ratios" % tp, fontsize=9, loc="left")
        ax.spines[["top", "right"]].set_visible(False)
        # (b) permutation null
        ax = AX[i, 1]
        ax.hist(r["null"], bins=40, color="#b0b7bd", edgecolor="white", linewidth=0.3)
        ax.axvline(r["aucR"], color=CR, lw=1.6)
        ax.annotate("observed\nAUC = %.3f\n$p$ = %.1e" % (r["aucR"], r["p"]),
                    xy=(r["aucR"], ax.get_ylim()[1] * 0.92), xytext=(-6, 0),
                    textcoords="offset points", ha="right", va="top", fontsize=7, color=CR)
        ax.set_xlabel("LOO-LDA AUC (66 ratios $\\rightarrow$ PCA15)")
        ax.set_ylabel("permutations")
        ax.set_title("%s  (b) label permutation (B=%d)" % (tp, B_PERM), fontsize=9, loc="left")
        ax.spines[["top", "right"]].set_visible(False)
        # (c) dose-response
        ax = AX[i, 2]
        S = r["S"]
        gs = [("negctrl", -0.18, "#7f8c8d"), ("sham", 0.0, CB), ("mix25", 0.25, "#e8a0a0"),
              ("mix50", 0.50, "#d96c6c"), ("mix75", 0.75, "#c44"), ("blast", 1.0, CR)]
        for g, x, c in gs:
            if g not in S: continue
            v = S[g]
            ax.scatter(np.full(len(v), x) + np.random.default_rng(SEED).normal(0, 0.012, len(v)),
                       v, s=9, color=c, alpha=0.55, edgecolor="none", zorder=2)
            ax.errorbar([x], [v.mean()], yerr=[sem(v)], fmt="o", ms=5, color=c,
                        mec="k", mew=0.5, capsize=3, lw=1.0, zorder=3)
        xs = [x for g, x, c in gs if g in S and g != "negctrl"]
        ms = [S[g].mean() for g, x, c in gs if g in S and g != "negctrl"]
        ax.plot(xs, ms, "-", color="k", lw=0.8, alpha=0.5, zorder=1)
        ax.axhline(0, c="k", lw=0.4, ls=":"); ax.axhline(1, c="k", lw=0.4, ls=":")
        ax.set_xticks([-0.18, 0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xticklabels(["neg\nctrl", "sham", "25%", "50%", "75%", "blast"], fontsize=7)
        ax.set_ylabel("peak-ratio axis score\n(sham=0, blast=1)")
        ax.set_title("%s  (c) dilution dose-response\n$\\rho_{mix}$=%+.2f  $\\rho_{all}$=%+.2f"
                     % (tp, r["rho_mix"], r["rho_all"]), fontsize=9, loc="left")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("exp2 chirality peak features: sham vs blast-conditioned media "
                 "(top row 2 h, bottom row 24 h)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    for ext in ("png", "svg"):
        p = os.path.join(HERE, "exp2_chirality.%s" % ext)
        fig.savefig(p, dpi=300 if ext == "png" else None)
        pr("saved " + p)

# ---------------------------------------------------------------- main
if __name__ == "__main__":
    pr("exp2 chirality peak-feature classification -- seed=%d, B_perm=%d" % (SEED, B_PERM))
    pr("caches: %s/exp2_{2h,24h}_peakfeats.npz" % CACHE)
    res = [run("2h"), run("24h")]
    pr("\n" + "#" * 78)
    pr("####  SUMMARY -- LOO-LDA AUC, sham vs blast")
    pr("#" * 78)
    pr("  %-22s %10s %10s %10s" % ("feature set", "exp2 2h", "exp2 24h", "exp3 3H"))
    pr("  %-22s %10.3f %10.3f %10s" % ("12 AUC fractions (PCA8)", res[0]["aucF"], res[1]["aucF"], "0.613"))
    pr("  %-22s %10.3f %10.3f %10s" % ("66 log-ratios (PCA15)", res[0]["aucR"], res[1]["aucR"], "0.938"))
    pr("  %-22s %10.3f %10.3f %10s" % ("78 combined (PCA15)", res[0]["aucB"], res[1]["aucB"], "0.920"))
    pr("  %-22s %10s %10s %10s" % ("full shape vector", "1.000", "1.000", "0.903"))
    pr("\n  66-ratio permutation p: 2h %.2e (k=%d) | 24h %.2e (k=%d)" % (
        res[0]["p"], res[0]["k"], res[1]["p"], res[1]["k"]))
    pr("  rho_mix: 2h %+.3f | 24h %+.3f   (shape-vector reference 0.86-0.94 at 24 h)" % (
        res[0]["rho_mix"], res[1]["rho_mix"]))
    figure(res)
    with open(os.path.join(HERE, "exp2_chirality_results.txt"), "w") as f:
        f.write("\n".join(LOG) + "\n")
    print("\nwrote", os.path.join(HERE, "exp2_chirality_results.txt"))
