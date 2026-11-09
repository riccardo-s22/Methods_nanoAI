"""Compare exp2 2h vs 24h sham-vs-TBI shape-difference spectra (same plate, two timepoints).
Is the 24h signature the SAME spectral change as 2h (matured), or a different/inverted one?
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
import numpy as np, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
import exp2_load
def pr(*a): print(*a); sys.stdout.flush()

DIRS = {"2h": "/mnt/c/Users/riccardo-s/Documents/CNT/TBI/exp2/2h",
        "24h": "/mnt/c/Users/riccardo-s/Documents/CNT/TBI/exp2/24h"}

def diff_spectrum(d):
    exp2_load.DIR = d
    cube, em = exp2_load.load_all()
    rows, cols, grp = exp2_load.labels()
    def shp(C):
        S = C.copy()
        for k in range(C.shape[2]): S[:, :, k] /= S[:, :, k].sum(1, keepdims=True)
        return S
    ds = shp(cube[grp == "TBI"]).mean(0) - shp(cube[grp == "sham"]).mean(0)  # (512,5)
    return em, ds.mean(1), ds  # mean over excitations + full

em, d2, d2f = diff_spectrum(DIRS["2h"])
_,  d24, d24f = diff_spectrum(DIRS["24h"])

band = (em >= 900) & (em <= 1450)
def cos(a, b): return float(a@b / (np.linalg.norm(a)*np.linalg.norm(b)))
pr("cosine(2h diff, 24h diff)  full mean-spectrum, 900-1450nm = %+.3f" % cos(d2[band], d24[band]))
pr("per-excitation cosines:")
for k, ex in enumerate(exp2_load.EXC):
    pr("  ex%3d  cos=%+.3f" % (ex, cos(d2f[band, k], d24f[band, k])))

# centroid of positive (gain) vs negative (loss) lobes
def lobe(em, d):
    g = d.copy(); g[g < 0] = 0; l = -d.copy(); l[l < 0] = 0
    return (em*g).sum()/g.sum(), (em*l).sum()/l.sum()
for tag, d in [("2h", d2), ("24h", d24)]:
    gc, lc = lobe(em[band], d[band])
    pr("%s: gain-lobe centroid %.0f nm,  loss-lobe centroid %.0f nm  (gain-loss=%+.0f nm)" % (tag, gc, lc, gc-lc))

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(em, d2, lw=2, c="#e67e22", label="2h  (TBI−sham)")
ax.plot(em, d24, lw=2, c="#cb4335", label="24h (TBI−sham)")
ax.axhline(0, c="grey", ls="--", lw=.8); ax.set_xlim(900, 1450)
ax.set_xlabel("emission (nm)"); ax.set_ylabel("TBI − sham (shape-normalized, mean over exc)")
ax.set_title("exp2 — sham→TBI shape-difference spectrum: 2h vs 24h (same plate)\n"
             "cos(2h,24h)=%+.2f" % cos(d2[band], d24[band]))
ax.legend()
fig.tight_layout()
OUT = os.environ.get("EXP2_OUT", "/mnt/c/Users/riccardo-s/Documents/CNT/targetALS/exp5/tbi_sham_tbi_outputs")
os.makedirs(OUT, exist_ok=True)
p = os.path.join(OUT, "exp2_2h_vs_24h_diffspectrum.png"); fig.savefig(p, dpi=150)
pr("saved", p)
