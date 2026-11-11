"""
exp3 loader — SWCNT exposed to brain-organoid conditioned media.

Format (per README): 2 plates x 2 timepoints (3H,6H) x 5 excitations (570,640,670,750,780).
Each .txt: 512 emission rows x 98 cols.
  col0 = emission axis (nm), col1 = synthetic ramp marker (excluded),
  cols 2..97 = 96 wells, ROW-MAJOR. well i -> row=i//12 (A-H), plate-col=i%12+1 (1-12).

Plate layout (README):
  SHAM plate: sham=cols1-4, resveratrol=cols5-8, MG132=cols9-12.
     exceptions (sham plate only): G9-12 & H9-10 = SWCNT in NON-conditioned media;
     H11-12 = media without SWCNTs (empty).  -> MG132-sham = cols9-12 rows A-F only.
  TBI plate: resveratrol=cols1-4, TBI=cols5-8, MG132=cols9-12.

Filenames (irregular): SHAM_EX_3_3H__{ex}.txt, SHAM_EX3_6H__{ex}.txt,
                       TBI_EX3_3H__{ex}.txt, TBI_EX3_6H__{ex}.txt
"""
import numpy as np, glob, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from repo_paths import data_path, env_path

DIR = env_path("EXP3_DIR", data_path("04_drug_response_prediction", "spectra"))
EXC = [570, 640, 670, 750, 780]
ROWS = "ABCDEFGH"

def wi(rc):
    return ROWS.index(rc[0]) * 12 + (int(rc[1:]) - 1)

def _files(plate, tp):
    """plate in {'SHAM','TBI'}, tp in {'3H','6H'} -> {excitation: path}."""
    out = {}
    for f in glob.glob(os.path.join(DIR, "*.txt")):
        b = os.path.basename(f)
        if "__" not in b:
            continue
        if not b.upper().startswith(plate):
            continue
        if tp not in b:
            continue
        ex = int(b.split("__")[-1].split(".")[0])
        out[ex] = f
    return out

def load_plate(plate, tp):
    """-> cube (96 wells, 512 emission, 5 excitation), em(512)."""
    files = _files(plate, tp)
    assert set(files) == set(EXC), "missing excitations for %s %s: got %s" % (plate, tp, sorted(files))
    em = None
    cube = np.zeros((96, 512, len(EXC)))
    for k, ex in enumerate(EXC):
        t = np.loadtxt(files[ex], delimiter="\t")
        if em is None:
            em = t[:, 0]
        cube[:, :, k] = t[:, 2:].T   # cols 2..97 -> 96 wells (col1 = ramp, skipped)
    return cube, em

# ---- well-group index sets (row-major 0-based) ----
def groups(plate):
    g = {}
    if plate == "SHAM":
        g["sham"]      = [r*12+c for r in range(8) for c in range(0, 4)]     # cols1-4
        g["resv_sham"] = [r*12+c for r in range(8) for c in range(4, 8)]     # cols5-8
        g["mg132_sham"]= [r*12+c for r in range(6) for c in range(8, 12)]    # cols9-12 rows A-F
        g["swcnt_noncm"]= [wi(x) for x in ("G9","G10","G11","G12","H9","H10")]
        g["no_swcnt"]  = [wi("H11"), wi("H12")]
    elif plate == "TBI":
        # NOTE (2026-06-27): user confirmed the TBI-plate resveratrol/MG132 columns were SWAPPED
        # vs the README. Corrected to ACTUAL layout: MG132=cols1-4, TBI=cols5-8, resveratrol=cols9-12.
        g["mg132_tbi"] = [r*12+c for r in range(8) for c in range(0, 4)]     # cols1-4  (actually MG132)
        g["tbi"]       = [r*12+c for r in range(8) for c in range(4, 8)]     # cols5-8
        g["resv_tbi"]  = [r*12+c for r in range(8) for c in range(8, 12)]    # cols9-12 (actually resveratrol)
    return {k: np.array(v) for k, v in g.items()}

def row_of(idx):
    return np.array([i // 12 for i in idx])
def col_of(idx):
    return np.array([i % 12 + 1 for i in idx])

if __name__ == "__main__":
    np.set_printoptions(suppress=True, precision=2, linewidth=160)
    for plate in ("SHAM", "TBI"):
        for tp in ("3H", "6H"):
            cube, em = load_plate(plate, tp)
            print("\n=== %s %s  cube %s  em %.0f-%.0f nm ===" % (plate, tp, cube.shape, em.min(), em.max()))
            G = groups(plate)
            for name, idx in G.items():
                tot = cube[idx, :, 0].sum(1)
                print("  %-13s n=%2d  bright@570 med=%.2e  peak med=%.0fnm"
                      % (name, len(idx), np.median(tot), np.median([em[np.argmax(cube[i,:,0])] for i in idx])))
