"""
exp3 7-day timepoint loader — plate A (and later B).

Format identical to the main exp3 files: each .txt = 512 emission rows x 98 cols.
  col0 = emission axis (nm), col1 = synthetic ramp marker (excluded),
  cols 2..97 = 96 wells ROW-MAJOR. well i -> row=i//12 (A-H), plate-col=i%12+1 (1-12).

Plate A column layout (user, 2026-07-16) — SAME-PLATE interleaved design
(sham & TBI on the same plate -> condition is NOT confounded with plate/batch):
  1  sham            2  tbi            3  sham+mg132     4  sham+res
  5  tbi+mg132       6  tbi+res        7  sham           8  tbi
  9  sham           10  sham+res      11  tbi+mg132     12  tbi+res

Files: 7d/A/TBI_d7__{ex}.txt , 7d/B/TBI_d7B__{ex}.txt
"""
import numpy as np, glob, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "shared"))
from repo_paths import data_path, env_path

BASE = env_path("EXP3_DIR", data_path("04_drug_response_prediction", "spectra"))
EXC = [570, 640, 670, 750, 780]
ROWS = "ABCDEFGH"

# plate-column (1-based) -> condition/drug arm
PLATE_A_COLS = {
    1: "sham",       2: "tbi",        3: "mg132_sham", 4: "resv_sham",
    5: "mg132_tbi",  6: "resv_tbi",   7: "sham",       8: "tbi",
    9: "sham",      10: "resv_sham", 11: "mg132_tbi", 12: "resv_tbi",
}
# Plate B layout (user, 2026-07-16) — clean balanced design, 16 wells/arm:
PLATE_B_COLS = {
    1: "sham",       2: "tbi",        3: "mg132_sham", 4: "resv_sham",
    5: "mg132_tbi",  6: "resv_tbi",   7: "sham",       8: "tbi",
    9: "mg132_sham",10: "resv_sham", 11: "mg132_tbi", 12: "resv_tbi",
}
LAYOUTS = {"A": PLATE_A_COLS, "B": PLATE_B_COLS}

def layout_for(plate):
    return LAYOUTS[plate]

def _files(plate):
    """plate in {'A','B'} -> {excitation: path}."""
    d = os.path.join(BASE, "7d", plate)
    out = {}
    for f in glob.glob(os.path.join(d, "*.txt")):
        b = os.path.basename(f)
        if "__" not in b:
            continue
        ex = int(b.split("__")[-1].split(".")[0])
        out[ex] = f
    return out

def load_plate(plate="A"):
    """-> cube (nwell, 512 emission, 5 excitation), em(512).
    nwell = ncols-2 (col0=emission, col1=ramp). Plate A has 96; plate B has 90
    (row-major acquisition, last 6 wells H7-H12 not captured -> wells A1..H6 only)."""
    files = _files(plate)
    assert set(files) == set(EXC), "missing excitations for 7d/%s: got %s" % (plate, sorted(files))
    em = None; cube = None
    for k, ex in enumerate(EXC):
        t = np.loadtxt(files[ex], delimiter="\t")
        if em is None:
            em = t[:, 0]; cube = np.zeros((t.shape[1] - 2, 512, len(EXC)))
        assert t.shape[1] - 2 == cube.shape[0], "well count mismatch across excitations for 7d/%s" % plate
        cube[:, :, k] = t[:, 2:].T   # cols 2.. -> wells row-major (col1 = ramp, skipped)
    return cube, em

def groups(layout=PLATE_A_COLS, nwell=96, keep=None):
    """condition/arm -> array of row-major 0-based well indices (all 8 rows per column),
    filtered to indices that actually exist (< nwell) so partial plates (B: H7-H12 missing) work.
    If keep (bool array over wells) is given, drop wells where keep is False (outlier removal)."""
    g = {}
    for col1, name in layout.items():
        idx = [r * 12 + (col1 - 1) for r in range(8) if r * 12 + (col1 - 1) < nwell]
        if keep is not None:
            idx = [i for i in idx if keep[i]]
        g.setdefault(name, []).extend(idx)
    return {k: np.array(sorted(v)) for k, v in g.items()}

def col_of(idx):
    return np.array([i % 12 + 1 for i in idx])
def row_of(idx):
    return np.array([i // 12 for i in idx])

def _robust_z(x):
    med = np.median(x); mad = np.median(np.abs(x - med))
    return (x - med) / (1.4826 * mad + 1e-30)

def outlier_keep(cube, dead_frac=0.15, z_intens=3.5, z_spec=5.0):
    """Boolean keep-mask over wells. Three-stage (see exp3_7d_outliers.py):
    (D) dead: total<=0 or <dead_frac*median; (I) position-detrended log-intensity |MAD-z|>z_intens;
    (S) shape-PCA(1..5) |MAD-z|>z_spec. Stages 2-3 fit on survivors of stage 1 so dead wells'
    garbage spectra don't distort them. Returns keep(bool, N)."""
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    N = cube.shape[0]
    row = np.array([i // 12 for i in range(N)], float); col = np.array([i % 12 for i in range(N)], float)
    tot = cube.sum(axis=(1, 2))
    dead = (tot <= 0) | (tot < dead_frac * np.median(tot))
    ok = ~dead
    Xd = np.column_stack([np.ones(ok.sum()), row[ok], col[ok]])
    beta, *_ = np.linalg.lstsq(Xd, np.log(tot[ok]), rcond=None)
    zi = _robust_z(np.log(tot[ok]) - Xd @ beta)
    intens = np.zeros(N, bool); intens[np.where(ok)[0]] = np.abs(zi) > z_intens
    X = cube.astype(float); X = X / np.clip(X.sum(1, keepdims=True), 1e-30, None); X = X.reshape(N, -1)
    Z = PCA(min(8, ok.sum() - 1), random_state=0).fit_transform(StandardScaler().fit_transform(X[ok]))
    zs = np.vstack([_robust_z(Z[:, k]) for k in range(min(5, Z.shape[1]))]).T
    spec = np.zeros(N, bool); spec[np.where(ok)[0]] = (np.abs(zs) > z_spec).any(1)
    return ~(dead | intens | spec)

if __name__ == "__main__":
    np.set_printoptions(suppress=True, precision=2, linewidth=160)
    cube, em = load_plate("A")
    print("=== 7d plate A  cube %s  em %.0f-%.0f nm ===" % (cube.shape, em.min(), em.max()))
    G = groups()
    for name, idx in G.items():
        tot = cube[idx, :, 0].sum(1)
        peak = np.median([em[np.argmax(cube[i, :, 0])] for i in idx])
        print("  %-12s n=%2d cols=%s  bright@570 med=%.2e  peak med=%.0fnm"
              % (name, len(idx), sorted(set(col_of(idx))), np.median(tot), peak))
