"""
exp2/2h loader. 5 files = 5 excitations (570,640,670,750,780).
Each file: 512 emission rows x 97 cols. col0=emission axis, cols1..96 = 96 wells row-major.
Well index i (0-based) -> plate row=i//12 (A-H), col=i%12+1 (1-12).
A1 (index0) = synthetic 0->1000 ramp placeholder -> excluded.

NEW interleaved layout (exp2):
  col1 = negctrl (SWCNT in media)
  cols2-9: even(2,4,6,8)=sham, odd(3,5,7,9)=TBI
  col10 = 3/4 TBI + 1/4 sham ; col11 = 1/2+1/2 ; col12 = 1/4 TBI + 3/4 sham
"""
import numpy as np, glob, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from repo_paths import data_path, env_path

# directory selectable via env EXP2_DIR (default 2h timepoint);
# loader is filename-agnostic (globs *.txt, parses excitation after "__").
# Default: <repo>/data/01_swcnt_blast_classification/spectra/2h, whatever EXP2_TAG says.
# EXP2_TAG only labels titles and output file names: for the 24 h read set BOTH
# EXP2_TAG=24h and EXP2_DIR=<repo>/data/01_swcnt_blast_classification/spectra/24h.
DIR = env_path("EXP2_DIR", data_path("01_swcnt_blast_classification", "spectra", "2h"))
EXC = [570, 640, 670, 750, 780]
ROWS = "ABCDEFGH"

def load_all():
    files = {int(f.split("__")[-1].split(".")[0]): f for f in glob.glob(os.path.join(DIR, "*.txt"))}
    em = None
    cube = np.zeros((96, 512, len(EXC)))  # well, emission, excitation
    for k, ex in enumerate(EXC):
        t = np.loadtxt(files[ex], delimiter="\t")
        if em is None: em = t[:, 0]
        cube[:, :, k] = t[:, 1:].T  # cols 1..96 -> wells
    return cube, em

def labels():
    rows, cols, grp = [], [], []
    for i in range(96):
        r, c = i // 12, i % 12 + 1
        rows.append(ROWS[r]); cols.append(c)
        if c == 1: g = "negctrl"
        elif c in (2, 4, 6, 8): g = "sham"
        elif c in (3, 5, 7, 9): g = "TBI"
        # NOTE (2026-06-16): user flagged cols 10/12 may be inverted.
        # Set MIX_INVERTED=1 in env to use col10=1/4 TBI, col12=3/4 TBI.
        elif c == 10: g = "mix25" if os.environ.get("MIX_INVERTED") else "mix75"
        elif c == 11: g = "mix50"
        elif c == 12: g = "mix75" if os.environ.get("MIX_INVERTED") else "mix25"
        grp.append(g)
    return np.array(rows), np.array(cols), np.array(grp)

if __name__ == "__main__":
    cube, em = load_all()
    rows, cols, grp = labels()
    # A1 ramp check
    print("A1 (index0) max=%.4g  -> ramp placeholder, excluded" % cube[0].max())
    valid = np.array([i for i in range(96) if i != 0])
    print("counts:", {g: int((grp[valid]==g).sum()) for g in ["negctrl","sham","TBI","mix75","mix50","mix25"]})
    # per-well total intensity (sum over emission, exc=570)
    tot = cube[:, :, 0].sum(1)
    print("\nper-well total @570 (8x12, A1=ramp):")
    np.set_printoptions(suppress=True, precision=2, linewidth=160)
    print(tot.reshape(8, 12) / 1e-11)
    # mean sham vs TBI total
    for ex_i, ex in enumerate(EXC):
        T = cube[:, :, ex_i].sum(1)
        s = T[(grp=="sham")]; t = T[(grp=="TBI")]
        print("ex%3d: sham tot %.3g  TBI tot %.3g  ratio %.3f" % (ex, s.mean(), t.mean(), t.mean()/s.mean()))
