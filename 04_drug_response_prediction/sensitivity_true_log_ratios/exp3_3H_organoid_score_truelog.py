"""
exp3 3H — per-organoid (= per-well) SWCNT injury score, for correlation with day-7 WB.

Axis = same as exp3_drug_axis.py: LDA sham-vs-TBI on 66 standardized log10 peak-amplitude ratios.
Honest scoring: sham & TBI wells are scored LEAVE-ONE-OUT (scaler+LDA refit without the well),
so their scores are not inflated by being in the training set. Drug-arm wells were never in
training -> scored with the full model. Axis oriented sham->TBI positive.

Residualized, WB-comparable score (like WB 'fc vs sham unt'):
    fc = 1 + (x - mean_sham) / (mean_TBI - mean_sham)      sham = 1, TBI = 2
(dividing by the sham mean itself is not usable: the LDA score is centred at ~0.)
Output: exp3_3H_organoid_scores.csv (one row per organoid/well).
"""
import os, sys, warnings
from itertools import combinations
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
warnings.filterwarnings("ignore")
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, ".."))                    # section top level: exp3_load.py
sys.path.insert(0, os.path.join(_HERE, "..", "..", "shared"))
from repo_paths import data_path, results_dir, env_path
import numpy as np, pandas as pd
import exp3_load as L
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.preprocessing import StandardScaler
OUT=env_path("EXP3_OUT", "") or results_dir("04_drug_response_prediction", "sensitivity_true_log_ratios"); os.makedirs(OUT, exist_ok=True)
CACHE_DIR=env_path("EXP3_CACHE_DIR", data_path("04_drug_response_prediction", "chirality_features"))  # shipped cache, read-only
CACHE=os.path.join(CACHE_DIR,"exp3_3H_peakfeats.npz")  # corrected labels
PAIRS=list(combinations(range(12),2)); flo=0.0  # TRUELOG RERUN (was 1e-6)
z=np.load(CACHE)
def ratio(g):
    M=z[f"MX_{g}"]; return np.column_stack([np.log10((M[:,a]+flo)/(M[:,b]+flo)) for a,b in PAIRS])
GR=["sham","tbi","resv_tbi","mg132_tbi","resv_sham","mg132_sham"]
PLATE={"sham":"SHAM","resv_sham":"SHAM","mg132_sham":"SHAM","tbi":"TBI","resv_tbi":"TBI","mg132_tbi":"TBI"}
COND={"sham":("sham","untreated"),"resv_sham":("sham","resveratrol"),"mg132_sham":("sham","MG132"),
      "tbi":("TBI","untreated"),"resv_tbi":("TBI","resveratrol"),"mg132_tbi":("TBI","MG132")}
R={g:ratio(g) for g in GR}
X=np.vstack([R["sham"],R["tbi"]]); y=np.r_[np.zeros(len(R["sham"])),np.ones(len(R["tbi"]))]

def fit(Xtr,ytr):
    sc=StandardScaler().fit(Xtr); Zt=sc.transform(Xtr); lda=LDA().fit(Zt,ytr)
    w=lda.coef_[0]; w=w*(1 if (Zt[ytr==1].mean(0)-Zt[ytr==0].mean(0))@w>0 else -1)
    wn=w/np.linalg.norm(w); return lambda A: sc.transform(A)@wn

full=fit(X,y)
score={g:full(R[g]) for g in GR}
loo=np.array([fit(np.delete(X,i,0),np.delete(y,i))(X[i:i+1])[0] for i in range(len(y))])
score["sham"],score["tbi"]=loo[y==0],loo[y==1]

ms,mt=score["sham"].mean(),score["tbi"].mean()
rows=[]
for g in GR:
    idx=L.groups(PLATE[g])[g]
    for k,(i,s) in enumerate(zip(idx,score[g])):
        rows.append(dict(plate=PLATE[g],well=L.ROWS[i//12]+str(i%12+1),injury=COND[g][0],treatment=COND[g][1],
                         group=g,axis_score=s,scoring="LOO" if g in("sham","tbi") else "full-model",
                         swcnt_fc=1+(s-ms)/(mt-ms)))
df=pd.DataFrame(rows)
p=os.path.join(OUT,"exp3_3H_organoid_scores.csv"); df.to_csv(p,index=False,float_format="%.5f")

from sklearn.metrics import roc_auc_score
print("reference means (LOO): sham axis=%.3f  TBI axis=%.3f   LOO AUC sham vs TBI=%.3f"%(ms,mt,roc_auc_score(y,loo)))
print("%-11s %3s  %8s %8s  %8s"%("group","n","fc mean","fc sd","fc SEM"))
for g in GR:
    f=df.loc[df.group==g,"swcnt_fc"]; print("%-11s %3d  %8.3f %8.3f  %8.3f"%(g,len(f),f.mean(),f.std(),f.sem()))
print("saved",p)
