"""Final supplementary-figure computation.
Experimental unit = one organoid = one well (n=32 sham, 32 blast).
Primary: leave-one-well-out. Controls: leave-one-plate-row-out, leave-one-column-out.
PCA is refit INSIDE every fold (no leakage)."""
import os, sys, json, itertools
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"): os.environ[v]="1"
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, ".."))                    # exp2_load.py
sys.path.insert(0, os.path.join(_HERE, "..", "..", "shared"))
from repo_paths import data_path, results_dir
import numpy as np, warnings; warnings.filterwarnings("ignore")
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.base import clone
from sklearn.model_selection import LeaveOneOut, LeaveOneGroupOut
from sklearn.metrics import roc_auc_score, roc_curve

P20 = lambda est: make_pipeline(StandardScaler(), PCA(20, random_state=0), est)
CLFS = {"LDA (PCA20)":P20(LDA()),
        "LogReg L2 (PCA20)":P20(LogisticRegression(max_iter=5000,C=1.0)),
        "Linear SVM (PCA20)":P20(SVC(kernel="linear",C=0.5)),
        "Elastic-net (PCA20)":P20(LogisticRegression(penalty="elasticnet",solver="saga",l1_ratio=.5,C=1.0,max_iter=20000)),
        "Elastic-net (full 2560)":make_pipeline(StandardScaler(),LogisticRegression(penalty="elasticnet",solver="saga",l1_ratio=.5,C=.5,max_iter=20000))}
NPERM=1000
out={}
for TAG,D in [("2h",data_path("01_swcnt_blast_classification","spectra","2h")),
              ("24h",data_path("01_swcnt_blast_classification","spectra","24h"))]:
    os.environ["EXP2_DIR"]=D
    for m in [k for k in list(sys.modules) if k.startswith("exp2_load")]: del sys.modules[m]
    from exp2_load import load_all, labels
    cube,em=load_all(); rows,cols,grp=labels()
    mask=np.isin(grp,["sham","TBI"])
    y=(grp[mask]=="TBI").astype(int); rr=rows[mask]; cc=cols[mask]; C=cube[mask]
    S=C.copy()
    for k in range(C.shape[2]): S[:,:,k]/=S[:,:,k].sum(1,keepdims=True)
    X=S.reshape(len(C),-1)
    r={"n_sham":int((y==0).sum()),"n_tbi":int((y==1).sum()),"y":y.tolist(),"auc":{},"roc":{},"scheme":{}}
    print(f"\n===== {TAG}: n={len(y)} organoids ({r['n_sham']} sham, {r['n_tbi']} blast)")
    def cv(model, splitter, groups=None):
        s=np.zeros(len(y))
        it = splitter.split(X,y,groups) if groups is not None else splitter.split(X)
        for tr,te in it: s[te]=clone(model).fit(X[tr],y[tr]).decision_function(X[te])
        return s
    # --- primary: leave-one-well-out, all classifiers
    for nm,mo in CLFS.items():
        s=cv(mo, LeaveOneOut())
        a=roc_auc_score(y,s); fpr,tpr,_=roc_curve(y,s)
        r["auc"][nm]=float(a); r["roc"][nm]=[fpr.tolist(),tpr.tolist()]
        if nm=="LDA (PCA20)": r["lda_loo_scores"]=s.tolist()
        print("   well-out  %-24s AUC=%.4f  gap=%+.1f" % (nm,a,s[y==1].min()-s[y==0].max()))
    # --- controls: row-out and column-out, all classifiers (PCA inside fold)
    for sch,gp,lab in [("row",rr,"plate row"),("col",cc,"plate column")]:
        d={}
        for nm,mo in CLFS.items():
            s=cv(mo, LeaveOneGroupOut(), gp)
            d[nm]={"auc":float(roc_auc_score(y,s)),"gap":float(s[y==1].min()-s[y==0].max())}
            print("   %s-out   %-24s AUC=%.4f  gap=%+.1f" % (sch,nm,d[nm]["auc"],d[nm]["gap"]))
        r["scheme"][sch]=d
    # --- permutation null, well level (primary), LDA
    s=cv(CLFS["LDA (PCA20)"], LeaveOneOut()); obs=roc_auc_score(y,s)
    Zpre=None
    rng=np.random.default_rng(0); null=[]
    for _ in range(NPERM):
        yp=rng.permutation(y)
        sp=np.zeros(len(yp))
        for tr,te in LeaveOneOut().split(X):
            sp[te]=clone(CLFS["LDA (PCA20)"]).fit(X[tr],yp[tr]).decision_function(X[te])
        null.append(roc_auc_score(yp,sp))
    null=np.array(null); c=int((null>=obs).sum())
    r["null_well"]=null.tolist(); r["perm_p_well"]=(c+1)/(NPERM+1); r["obs"]=float(obs)
    print("   PERMUTATION (well level, n=%d): obs=%.4f  null mean %.3f max %.3f  %d>=obs  p=%.4f"
          % (NPERM,obs,null.mean(),null.max(),c,(c+1)/(NPERM+1)))
    out[TAG]=r
JSON=os.path.join(results_dir("01_swcnt_blast_classification","supp_fig2_auc"),"supp_auc_results.json")
json.dump(out,open(JSON,"w")); print("\nwrote",JSON)
