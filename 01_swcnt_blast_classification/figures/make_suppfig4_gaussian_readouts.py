"""Supplementary Figure 4, Gaussian-fit version: same two readouts of the dilution series (exp2, 24 h read),
computed on the split-Gaussian chirality features instead of the shape-normalized spectra.

Features: 66 pairwise log10 ratios of the fitted split-Gaussian maxima (MX) of the 12 chiralities
(exactly as in paper/svg/exp2_chirality/exp2_chirality_classify.py: sham/blast MX from the cached
targetALS/exp5/exp2_24h_peakfeats.npz, dilution + negative-control MX from exp2_24h_mix_peakfeats.npz,
MIX_INVERTED=1 labelling).
  a  classifier score: StandardScaler -> PCA(15) -> LDA fitted on sham + blast (the chirality model),
     rescaled sham mean = 0, blast mean = 1
  b  projection onto d = mu_blast - mu_sham in the raw 66-ratio space (no standardization, no PCA)
  c  schematic
  d  cosine between each group's mean change from sham and d, same 66 features
Run:  python make_suppfig4_gaussian_readouts.py [out_basename]
"""
import os, sys
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[v] = "1"
from itertools import combinations
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.pipeline import make_pipeline
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Arc
plt.rcParams.update({"svg.fonttype": "none", "font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 8,
                     "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False})

# ------------------------------------------------------------------ data (Gaussian-fit chirality features)
CACHE = "/mnt/c/Users/riccardo-s/Documents/CNT/targetALS/exp5/exp2_24h_peakfeats.npz"
MIXC = "/mnt/c/Users/riccardo-s/Documents/CNT/TBI/paper/svg/exp2_chirality/exp2_24h_mix_peakfeats.npz"
PAIRS = list(combinations(range(12), 2))
# FEAT=paper  : log10((MXa + 1e-6) / (MXb + 1e-6)), exactly as in exp2_chirality_classify.py. NOTE: the fitted
#               maxima are ~1e-14..1e-12, so the 1e-6 floor dominates and these features are, to r = 1.0000,
#               linear differences of absolute peak heights (MXa - MXb), not intensity ratios.
# FEAT=truelog: log10(MXa / MXb), no floor (all MX > 0 in these caches).
FEAT = os.environ.get("FEAT", "paper")
flo = 1e-6 if FEAT == "paper" else 0.0
def ratios(MX):
    return np.column_stack([np.log10((MX[:, a] + flo) / (MX[:, b] + flo)) for a, b in PAIRS])
z = np.load(CACHE); zm = np.load(MIXC)
MXg = {"sham": z["MXs"], "TBI": z["MXt"], "mix25": zm["mix25"], "mix50": zm["mix50"],
       "mix75": zm["mix75"], "negctrl": zm["negctrl"]}
ORDER = ["sham", "mix25", "mix50", "mix75", "TBI", "negctrl"]
LAB = {"sham": "Sham", "mix25": "25% blast", "mix50": "50% blast", "mix75": "75% blast",
       "TBI": "Blast", "negctrl": "Non-conditioned\nmedium"}
COL = {"sham": "#2e86c1", "mix25": "#a9dfbf", "mix50": "#52be80", "mix75": "#1e8449",
       "TBI": "#cb4335", "negctrl": "#7f8c8d"}
MRK = {g: ("D" if g == "negctrl" else "o") for g in ORDER}
X = {g: ratios(MXg[g]) for g in ORDER}
print("FEAT =", FEAT); print("wells:", {g: len(X[g]) for g in ORDER})
Xtr = np.vstack([X["sham"], X["TBI"]]); ytr = np.r_[np.zeros(len(X["sham"])), np.ones(len(X["TBI"]))]

pipe = make_pipeline(StandardScaler(), PCA(n_components=15, random_state=0), LDA()).fit(Xtr, ytr)
raw = {g: pipe.decision_function(X[g]) for g in ORDER}
s0, s1 = raw["sham"].mean(), raw["TBI"].mean()
score = {g: (raw[g] - s0) / (s1 - s0) for g in ORDER}

mu_s, mu_b = X["sham"].mean(0), X["TBI"].mean(0); d = mu_b - mu_s; dn = d / (d @ d)
proj = {g: (X[g] - mu_s) @ dn for g in ORDER}
cos = {g: float(((X[g].mean(0) - mu_s) @ d) / (np.linalg.norm(X[g].mean(0) - mu_s) * np.linalg.norm(d)))
       for g in ORDER if g != "sham"}
print("rescaled LDA score:", {g: round(float(score[g].mean()), 2) for g in ORDER})
print("projection on d   :", {g: round(float(proj[g].mean()), 2) for g in ORDER})
print("cosine to d       :", {g: round(c, 2) for g, c in cos.items()})

# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(7.2, 6.4))                     # 183 mm wide (full page)
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], hspace=0.55, wspace=0.30,
                      left=0.14, right=0.93, top=0.88, bottom=0.09)
axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
axC, axD = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
ypos = {"sham": 0, "mix25": 1, "mix50": 2, "mix75": 3, "TBI": 4, "negctrl": 5.4}

def strips(ax, vals, xlabel, title, subtitle, show_labels, xlim=(-0.6, 1.8)):
    rng = np.random.default_rng(0)
    for xx, lab in [(0, "sham mean"), (1, "blast mean")]:
        ax.axvline(xx, color="0.55", lw=0.8, ls=(0, (3, 2)), zorder=1)
    for g in ORDER:
        v = vals[g]; y = ypos[g] + rng.uniform(-0.22, 0.22, len(v))
        ax.scatter(v, y, s=14, marker=MRK[g], c=COL[g], edgecolor="0.15", lw=0.35, alpha=0.9, zorder=3)
        m = float(v.mean())
        ax.plot([m, m], [ypos[g] - 0.38, ypos[g] + 0.38], color="k", lw=1.6, zorder=4, solid_capstyle="butt")
        ax.text(1.02, ypos[g], "%.2f" % (0.0 if abs(m) < 0.005 else m), transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=7.5)
    ax.text(1.02, -0.95, "mean", transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=7, color="0.35")
    ax.set_yticks([ypos[g] for g in ORDER]); ax.set_yticklabels([LAB[g] for g in ORDER] if show_labels else [])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(6.0, -1.3); ax.set_xlim(*xlim)
    if xlim == (-0.6, 1.8): ax.set_xticks([-0.5, 0, 0.5, 1, 1.5])
    ax.set_xlabel(xlabel)
    ax.text(0, -1.25, "sham\nmean" if xlim[1] - xlim[0] < 4 else "sham", ha="center", va="bottom", fontsize=6.5, color="0.35")
    ax.text(1, -1.25, "blast\nmean" if xlim[1] - xlim[0] < 4 else "blast", ha="center", va="bottom", fontsize=6.5, color="0.35")
    ax.text(-0.02, 1.25, title, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left")
    ax.text(-0.02, 1.17, subtitle, transform=ax.transAxes, fontsize=7.2, color="0.25", va="bottom", ha="left")

def lim(v):
    lo = min(-0.6, min(float(a.min()) for a in v.values()) - 0.2); hi = max(1.8, max(float(a.max()) for a in v.values()) + 0.2)
    return (np.floor(lo * 2) / 2, np.ceil(hi * 2) / 2) if (lo, hi) != (-0.6, 1.8) else (lo, hi)
strips(axA, score, "Classifier score\n(0 = sham mean, 1 = blast mean)",
       "a   Trained classifier score", "LDA on 66 Gaussian peak ratios, trained on sham vs blast only", True, lim(score))
strips(axB, proj, "Position along the sham→blast peak-ratio difference\n(0 = sham mean, 1 = blast mean)",
       "b   Position along sham→blast direction", "No model: projection of each well onto d (see c)", False, lim(proj))

# c: schematic of the projection and of the direction match
axC.set_xlim(-0.15, 1.25); axC.set_ylim(-0.30, 0.95); axC.set_aspect("equal"); axC.axis("off")
S0, B1 = np.array([0.0, 0.0]), np.array([1.0, 0.0]); P = np.array([0.62, 0.42])
axC.add_patch(FancyArrowPatch(S0, B1, arrowstyle="-|>", mutation_scale=11, color="k", lw=1.6, zorder=3))
axC.add_patch(FancyArrowPatch(S0, P, arrowstyle="-|>", mutation_scale=10, color=COL["mix50"], lw=1.6, zorder=3))
axC.plot([P[0], P[0]], [P[1], 0], color="0.4", lw=0.9, ls=(0, (2, 2)), zorder=2)
axC.plot([0, P[0]], [-0.10, -0.10], color="0.25", lw=0.9); axC.plot([0, 0], [-0.13, -0.07], color="0.25", lw=0.9)
axC.plot([P[0], P[0]], [-0.13, -0.07], color="0.25", lw=0.9)
axC.text(P[0] / 2, -0.15, "projection (panel b)", ha="center", va="top", fontsize=7)
ang = np.degrees(np.arctan2(P[1], P[0]))
axC.add_patch(Arc(S0, 0.42, 0.42, theta1=0, theta2=ang, color="0.25", lw=0.9))
axC.text(0.26, 0.055, "θ", fontsize=8)
for pt, g, lab, dx, ha in [(S0, "sham", "sham mean", -0.02, "right"), (B1, "TBI", "blast mean", 0.02, "left")]:
    axC.scatter(*pt, s=40, c=COL[g], edgecolor="k", lw=0.6, zorder=4)
    axC.text(pt[0] + dx, pt[1] + 0.06, lab, ha=ha, fontsize=7)
axC.scatter(*P, s=40, c=COL["mix50"], edgecolor="k", lw=0.6, zorder=4)
axC.text(P[0] + 0.03, P[1] + 0.02, "peak ratios of a sample", fontsize=7, va="bottom")
axC.text(1.0, -0.035, "d", fontsize=8, fontstyle="italic", ha="right", va="top")
axC.text(-0.12, 0.95, "d = mean blast − mean sham peak ratios\nprojection = distance moved along d\n   (sham = 0, blast = 1)\ndirection match = cos θ\n   (1 = change points toward blast)", fontsize=6.8, va="top", color="0.25")
axC.set_title("c   How panels b and d are computed", loc="left", fontsize=9, fontweight="bold", pad=8)

# d: direction match of each group's mean change with d
gD = ["mix25", "mix50", "mix75", "TBI", "negctrl"]; yD = {g: i for i, g in enumerate(gD)}
for g in gD:
    axD.plot([0, cos[g]], [yD[g], yD[g]], color="0.75", lw=1.2, zorder=1)
    axD.scatter(cos[g], yD[g], s=46, marker=MRK[g], c=COL[g], edgecolor="k", lw=0.6, zorder=3)
    axD.text(cos[g] + (0.04 if cos[g] >= 0 else -0.04), yD[g], "%.2f" % cos[g], va="center", ha="left" if cos[g] >= 0 else "right", fontsize=7.5)
axD.axvline(0, color="0.55", lw=0.8)
axD.set_yticks([yD[g] for g in gD]); axD.set_yticklabels([LAB[g] for g in gD]); axD.tick_params(axis="y", length=0)
axD.set_ylim(len(gD) - 0.4, -0.6); cmin = min(cos.values()); axD.set_xlim(min(-0.1, cmin - 0.35), 1.15)
if cmin >= 0: axD.set_xticks([0, 0.25, 0.5, 0.75, 1])
axD.set_xlabel("Direction match with blast\n(cos θ of the group-mean change from sham)")
axD.set_title("d   Direction of change versus blast", loc="left", fontsize=9, fontweight="bold", pad=8)

out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "suppfig4_gaussian_readouts" + ("" if FEAT == "paper" else "_truelog"))
fig.savefig(out + ".svg"); fig.savefig(out + ".png", dpi=300)
from PIL import Image
Image.open(out + ".png").convert("RGB").save(out + ".tif", dpi=(300, 300), compression="tiff_lzw")
print("saved", out + ".svg/.png/.tif")
