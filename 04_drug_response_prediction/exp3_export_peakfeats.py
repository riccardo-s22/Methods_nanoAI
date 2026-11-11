"""Export the canonical fitted-peak chirality features to CSV (both 3H & 6H) into the exp3 folder.
Reads the per-timepoint npz caches (exp3_{TP}_peakfeats.npz) written by the analysis scripts.
For each well: group, plate-row/col, 12 chirality gauss_max, 12 gauss_auc, 12 AUC-fractions,
and 66 pairwise log10 intensity ratios. One CSV per timepoint + a combined CSV."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from repo_paths import data_path, results_dir, env_path
from itertools import combinations
import exp3_load as L
OUT=env_path("EXP3_OUT", "") or results_dir("04_drug_response_prediction"); os.makedirs(OUT, exist_ok=True)
CACHE_DIR=env_path("EXP3_CACHE_DIR", data_path("04_drug_response_prediction", "chirality_features"))  # shipped cache, read-only
ORDER=["ch8_3","ch6_5","ch7_5","ch10_2","ch9_4","ch8_4","ch7_6","ch8_6","ch8_7","ch9_5","ch10_3","ch10_5"]
SH={k:k.replace("ch","").replace("_","_") for k in ORDER}   # 8_3, 6_5, ...
PAIRS=list(combinations(range(12),2)); PNAME=[f"{SH[ORDER[a]]}/{SH[ORDER[b]]}" for a,b in PAIRS]
GROUPS=["sham","tbi","resv_tbi","mg132_tbi","resv_sham","mg132_sham","swcnt_noncm"]
# plate + well-index for each group (to record row/col)
PLATE={"sham":"SHAM","resv_sham":"SHAM","mg132_sham":"SHAM","swcnt_noncm":"SHAM",
       "tbi":"TBI","resv_tbi":"TBI","mg132_tbi":"TBI"}
def gidx(g):
    Gs=L.groups("SHAM"); Gt=L.groups("TBI")
    return (Gs if PLATE[g]=="SHAM" else Gt)[ {"sham":"sham","resv_sham":"resv_sham","mg132_sham":"mg132_sham",
        "swcnt_noncm":"swcnt_noncm","tbi":"tbi","resv_tbi":"resv_tbi","mg132_tbi":"mg132_tbi"}[g] ]
flo=1e-6
def export(tp):
    cache=os.path.join(CACHE_DIR,"exp3_%s_peakfeats.npz"%tp)
    if not os.path.exists(cache): print("missing",cache); return None
    z=np.load(cache); rows=[]
    for g in GROUPS:
        MX=z[f"MX_{g}"]; AU=z[f"AU_{g}"]; idx=gidx(g)
        fr=AU/np.clip(AU.sum(1,keepdims=True),1e-30,None)
        for i in range(len(MX)):
            w=idx[i]; rec={"timepoint":tp,"group":g,"plate":PLATE[g],
                           "well_row":"ABCDEFGH"[w//12],"well_col":w%12+1}
            for k,nm in enumerate(ORDER):
                rec[f"max_{SH[nm]}"]=MX[i,k]; rec[f"auc_{SH[nm]}"]=AU[i,k]; rec[f"frac_{SH[nm]}"]=fr[i,k]
            for j,(a,b) in enumerate(PAIRS):
                rec[f"logratio_{PNAME[j]}"]=np.log10((MX[i,a]+flo)/(MX[i,b]+flo))
            rows.append(rec)
    df=pd.DataFrame(rows); p=os.path.join(OUT,"exp3_%s_chirality_peak_features.csv"%tp)
    df.to_csv(p,index=False); print("saved %s  (%d wells x %d cols)"%(p,len(df),df.shape[1])); return df
dfs=[d for d in (export("3H"),export("6H")) if d is not None]
if len(dfs)==2:
    allp=os.path.join(OUT,"exp3_chirality_peak_features_3H_6H.csv")
    pd.concat(dfs,ignore_index=True).to_csv(allp,index=False); print("saved combined",allp)
