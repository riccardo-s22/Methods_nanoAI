"""
exp3 3H — test the possible TBI-plate label swap (resveratrol cols1-4 <-> MG132 cols9-12).
The sham->TBI axis is unchanged (built from sham & TBI, neither swapped). Only the TBI-plate drug
arm labels swap. Report key drug-axis metrics under BOTH labelings side by side so the correct map
can be chosen from lab records (NOT from which fits the hypothesis).
  cache arrays: 'resv_tbi' = TBI cols1-4 ; 'mg132_tbi' = TBI cols9-12 (as loaded per README).
  ORIGINAL: resveratrol=cols1-4, MG132=cols9-12.   SWAPPED: resveratrol=cols9-12, MG132=cols1-4.
Sham plate fixed: resv_sham=cols5-8, mg132_sham=cols9-12.
"""
import os, sys, warnings
from itertools import combinations
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
warnings.filterwarnings("ignore")
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, ".."))                    # section top level: exp3_load.py
sys.path.insert(0, os.path.join(_HERE, "..", "..", "shared"))
from repo_paths import data_path, env_path
import numpy as np
from scipy import stats
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.preprocessing import StandardScaler
def pr(*a): print(*a); sys.stdout.flush()
ORDER=["ch8_3","ch6_5","ch7_5","ch10_2","ch9_4","ch8_4","ch7_6","ch8_6","ch8_7","ch9_5","ch10_3","ch10_5"]
PAIRS=list(combinations(range(12),2)); flo=1e-6
CACHE_DIR=env_path("EXP3_CACHE_DIR", data_path("04_drug_response_prediction", "chirality_features"))  # shipped cache, read-only
z=np.load(os.path.join(CACHE_DIR,"exp3_3H_peakfeats.npz"))
def ratio(g):
    M=z[f"MX_{g}"]; return np.column_stack([np.log10((M[:,a]+flo)/(M[:,b]+flo)) for a,b in PAIRS])
base={g:ratio(g) for g in ["sham","tbi","resv_tbi","mg132_tbi","resv_sham","mg132_sham"]}

# axis from sham vs tbi (unchanged by swap)
sc=StandardScaler().fit(np.vstack([base["sham"],base["tbi"]]))
Z={g:sc.transform(base[g]) for g in base}
lda=LDA().fit(np.vstack([Z["sham"],Z["tbi"]]),np.r_[np.zeros(len(Z["sham"])),np.ones(len(Z["tbi"]))])
w=lda.coef_[0]; w=(1 if (Z["tbi"].mean(0)-Z["sham"].mean(0))@w>0 else -1)*w; wn=w/np.linalg.norm(w)
score={g:Z[g]@wn for g in base}

def perm_t(a,b,B=10000,seed=0):
    obs=b.mean()-a.mean(); pool=np.r_[a,b]; n=len(a); rng=np.random.default_rng(seed); c=0
    for _ in range(B):
        p=rng.permutation(len(pool));
        if abs(pool[p][n:].mean()-pool[p][:n].mean())>=abs(obs): c+=1
    return obs,(c+1)/(B+1)
def did(armT,armS):
    a,b,c,d=score["tbi"],score[armT],score["sham"],score[armS]
    D=(b.mean()-a.mean())-(d.mean()-c.mean())
    se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b)+c.var(ddof=1)/len(c)+d.var(ddof=1)/len(d))
    return D,se,2*stats.t.sf(abs(D/se),len(a)+len(b)+len(c)+len(d)-4)
def pct(armT,seed=0):
    rng=np.random.default_rng(seed); b=[]
    for _ in range(5000):
        s=score["sham"][rng.integers(0,32,32)].mean(); t=score["tbi"][rng.integers(0,32,32)].mean()
        r=score[armT][rng.integers(0,len(score[armT]),len(score[armT]))].mean()
        if abs(t-s)>1e-9: b.append(100*(t-r)/(t-s))
    return np.median(b),np.percentile(b,2.5),np.percentile(b,97.5)

pr("axis means: sham=%.3f tbi=%.3f | TBI-cols1-4(=resv_tbi cache)=%.3f  TBI-cols9-12(=mg132_tbi cache)=%.3f"
   %(score["sham"].mean(),score["tbi"].mean(),score["resv_tbi"].mean(),score["mg132_tbi"].mean()))
pr("            resv_sham(cols5-8)=%.3f  mg132_sham(cols9-12)=%.3f"%(score["resv_sham"].mean(),score["mg132_sham"].mean()))

# labelings: name -> (resveratrol_TBI_arm, MG132_TBI_arm)
LAB={"ORIGINAL (README)":("resv_tbi","mg132_tbi"),   # resv=cols1-4, mg=cols9-12
     "SWAPPED":("mg132_tbi","resv_tbi")}              # resv=cols9-12, mg=cols1-4
for name,(rT,mT) in LAB.items():
    pr("\n================ %s : resveratrol=%s, MG132=%s ================"%(name,
        "TBI cols1-4" if rT=="resv_tbi" else "TBI cols9-12","TBI cols9-12" if mT=="mg132_tbi" else "TBI cols1-4"))
    for label,armT,armS in [("RESVERATROL",rT,"resv_sham"),("MG132",mT,"mg132_sham")]:
        o,p=perm_t(score["tbi"],score[armT])         # raw within-TBI-plate axis move
        D,se,pd=did(armT,armS)                        # sham-controlled DiD
        m,lo,hi=pct(armT)
        pr("  %-11s: within-TBI axis move=%+.3f (p=%.4f) | %% toward sham=%.0f%% (CI %.0f..%.0f) | DiD(sham-ctrl)=%+.3f p=%.4f"
           %(label,o,p,m,lo,hi,D,pd))
pr("\nNOTE: 'within-TBI axis move' has NO sham control -> flips with the swap; the DiD subtracts the")
pr("drug's own sham-plate effect. Choose the labeling from the PLATE MAP/lab records, not from fit.")
