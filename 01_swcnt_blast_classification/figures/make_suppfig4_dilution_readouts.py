"""Supplementary Figure 4 (rebuilt): two readouts of the dilution series, exp2 24 h read.

Same data and computations as make_concept_fig_realdata.py (classifier = StandardScaler -> PCA(20) -> LDA
fitted on sham + blast wells, as in Fig. 2D; model-free projection onto d = mu_blast - mu_sham in raw shape
space; cosine between each group's mean change from sham and d). Only the presentation changes:
  a  per-well classifier score, one row per group, on a 0 (sham mean) .. 1 (blast mean) scale
  b  per-well position along d (no model), same scale
  c  schematic of how b and d are computed
  d  direction match (cosine) of each group's mean spectral change with d
Writes <out>.svg (editable text), <out>.png and <out>.tif (300 dpi).
Run:  python make_suppfig4_dilution_readouts.py [out_basename]
"""
import os, sys
os.environ.setdefault("EXP2_DIR", "/mnt/c/Users/riccardo-s/Documents/CNT/TBI/exp2/24h")
os.environ.setdefault("MIX_INVERTED", "1")
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[v] = "1"
sys.path.insert(0, "/mnt/c/Users/riccardo-s/Documents/CNT/TBI/exp2")
import numpy as np
from exp2_load import load_all, labels
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Arc
plt.rcParams.update({"svg.fonttype": "none", "font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 8,
                     "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False})

# ------------------------------------------------------------------ data (unchanged from the original script)
cube, em = load_all(); rows, cols, grp = labels()
def shp(C):
    S = C.copy()
    for k in range(C.shape[2]): S[:, :, k] /= S[:, :, k].sum(1, keepdims=True)
    return S.reshape(len(C), -1)
def wells(g):
    idx = np.where(grp == g)[0]
    return idx[idx != 0] if g == "negctrl" else idx          # A1 is an export placeholder
ORDER = ["sham", "mix25", "mix50", "mix75", "TBI", "negctrl"]
LAB = {"sham": "Sham", "mix25": "25% blast", "mix50": "50% blast", "mix75": "75% blast",
       "TBI": "Blast", "negctrl": "Non-conditioned\nmedium"}
COL = {"sham": "#2e86c1", "mix25": "#a9dfbf", "mix50": "#52be80", "mix75": "#1e8449",
       "TBI": "#cb4335", "negctrl": "#7f8c8d"}
MRK = {g: ("D" if g == "negctrl" else "o") for g in ORDER}
X = {g: shp(cube[wells(g)]) for g in ORDER}
Xtr = np.vstack([X["sham"], X["TBI"]]); ytr = np.r_[np.zeros(len(X["sham"])), np.ones(len(X["TBI"]))]

sc = StandardScaler().fit(Xtr); pca = PCA(20).fit(sc.transform(Xtr))
Z = {g: pca.transform(sc.transform(X[g])) for g in ORDER}
lda = LDA().fit(np.vstack([Z["sham"], Z["TBI"]]), ytr)
raw = {g: lda.decision_function(Z[g]) for g in ORDER}
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

def strips(ax, vals, xlabel, title, subtitle, show_labels):
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
    ax.set_ylim(6.0, -1.3); ax.set_xlim(-0.6, 1.8); ax.set_xticks([-0.5, 0, 0.5, 1, 1.5])
    ax.set_xlabel(xlabel)
    ax.text(0, -1.25, "sham\nmean", ha="center", va="bottom", fontsize=6.5, color="0.35")
    ax.text(1, -1.25, "blast\nmean", ha="center", va="bottom", fontsize=6.5, color="0.35")
    ax.text(-0.02, 1.25, title, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left")
    ax.text(-0.02, 1.17, subtitle, transform=ax.transAxes, fontsize=7.2, color="0.25", va="bottom", ha="left")

strips(axA, score, "Classifier score\n(0 = sham mean, 1 = blast mean)",
       "a   Trained classifier score", "LDA trained on sham vs blast wells only (as in Figure 2D)", True)
strips(axB, proj, "Position along the sham→blast spectral difference\n(0 = sham mean, 1 = blast mean)",
       "b   Position along sham→blast direction", "No model: projection of each spectrum onto d (see c)", False)

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
axC.text(P[0] + 0.03, P[1] + 0.02, "spectrum of a sample", fontsize=7, va="bottom")
axC.text(1.0, -0.035, "d", fontsize=8, fontstyle="italic", ha="right", va="top")
axC.text(-0.12, 0.95, "d = mean blast − mean sham spectrum\nprojection = distance moved along d\n   (sham = 0, blast = 1)\ndirection match = cos θ\n   (1 = change points toward blast)", fontsize=6.8, va="top", color="0.25")
axC.set_title("c   How panels b and d are computed", loc="left", fontsize=9, fontweight="bold", pad=8)

# d: direction match of each group's mean change with d
gD = ["mix25", "mix50", "mix75", "TBI", "negctrl"]; yD = {g: i for i, g in enumerate(gD)}
for g in gD:
    axD.plot([0, cos[g]], [yD[g], yD[g]], color="0.75", lw=1.2, zorder=1)
    axD.scatter(cos[g], yD[g], s=46, marker=MRK[g], c=COL[g], edgecolor="k", lw=0.6, zorder=3)
    axD.text(cos[g] + 0.04, yD[g], "%.2f" % cos[g], va="center", fontsize=7.5)
axD.axvline(0, color="0.55", lw=0.8)
axD.set_yticks([yD[g] for g in gD]); axD.set_yticklabels([LAB[g] for g in gD]); axD.tick_params(axis="y", length=0)
axD.set_ylim(len(gD) - 0.4, -0.6); axD.set_xlim(-0.1, 1.15); axD.set_xticks([0, 0.25, 0.5, 0.75, 1])
axD.set_xlabel("Direction match with the blast spectrum\n(cos θ of the group-mean change from sham)")
axD.set_title("d   Direction of change versus blast", loc="left", fontsize=9, fontweight="bold", pad=8)

out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "suppfig4_dilution_readouts")
fig.savefig(out + ".svg"); fig.savefig(out + ".png", dpi=300)
from PIL import Image
Image.open(out + ".png").convert("RGB").save(out + ".tif", dpi=(300, 300), compression="tiff_lzw")
print("saved", out + ".svg/.png/.tif")
