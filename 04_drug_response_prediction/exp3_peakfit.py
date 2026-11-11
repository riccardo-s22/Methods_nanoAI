"""
exp3 — 3H ONLY, CANONICAL fitted-peak (chirality) features.
Reuses physical_descriptors/chirality_interpolated_extraction.py (the project's canonical
split-Gaussian-per-(n,m) descriptor extractor; same one used by tbi_chirality_gaussian.py),
with the 12 canonical CHIRALITY_POSITIONS. exp3 has only 5 excitations, so per well we build the
2D interpolated EEM surface over those 5 excitations and fit each chirality at its RESONANT
excitation slice (this is what separates the overlapping ~1130nm emission peaks 8_4/9_4/7_6).

Per well -> 12 chirality split-Gaussian fits (gauss_auc, gauss_max, gauss_em_center, gauss_fwhm).
Feature = per-well chirality AUC FRACTIONS (brightness-invariant chirality redistribution).
3H analysis: sham vs TBI (LOO-LDA AUC + permutation), which (n,m) drives it (Cohen d + t + BH-FDR),
center shifts, then project drug arms (resveratrol/MG132) onto the chirality-feature sham->TBI axis.
"""
import os, sys, warnings, importlib.util
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from repo_paths import data_path, results_dir, env_path
import numpy as np
import exp3_load as L
from scipy.stats import ttest_ind
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def pr(*a): print(*a); sys.stdout.flush()
OUT=env_path("EXP3_OUT", "") or results_dir("04_drug_response_prediction"); os.makedirs(OUT, exist_ok=True)

# ---- import canonical extractor ----
PD=os.path.join(os.path.dirname(os.path.abspath(__file__)),"vendored")
spec=importlib.util.spec_from_file_location("cie", os.path.join(PD,"chirality_interpolated_extraction.py"))
cie=importlib.util.module_from_spec(spec); spec.loader.exec_module(cie)
CHIR=cie.CHIRALITY_POSITIONS; ORDER=list(CHIR.keys())
SHORT={k:k.replace("ch","(").replace("_",",")+")" for k in ORDER}   # ch6_5 -> (6,5)
EXC=np.array([570,640,670,750,780],float)
COL={"sham":"#2e86c1","tbi":"#cb4335","resv_tbi":"#27ae60","mg132_tbi":"#8e44ad",
     "resv_sham":"#a9cce3","mg132_sham":"#d2b4de","swcnt_noncm":"#7f8c8d"}

# ---- load 3H plates ----
cs,em=L.load_plate("SHAM","3H"); ct,_=L.load_plate("TBI","3H")
Gs=L.groups("SHAM"); Gt=L.groups("TBI")
exmin,exmax=EXC.min(),EXC.max()

def well_chirality(mat):
    """mat (512 em, 5 ex) -> dict chir-> descriptors, mirroring cie.extract_sample_tp."""
    spline,em_f,ex_f=cie.build_interpolated_surface(mat, em, EXC, factor=cie.INTERP_FACTOR)
    out={}
    for name,c in CHIR.items():
        tex=min(max(c["excitation_nm"],exmin),exmax)         # clip to available excitations
        be,bx,mx=cie.find_local_max_in_radius(spline,em_f,ex_f,c["emission_nm"],tex,cie.SEARCH_RADIUS_NM)
        ep,ip=cie.extract_consistent_gradient_points(spline,em_f,be,bx)
        d=cie.fit_gaussian_and_extract(ep,ip,be,mx)
        out[name]=d
    return out

def group_feats(cube,idx):
    """-> AUC (n,12), CENTER (n,12), FWHM (n,12)."""
    A=np.zeros((len(idx),12)); C=np.zeros((len(idx),12)); Fw=np.zeros((len(idx),12)); R2=np.zeros((len(idx),12))
    for i,w in enumerate(idx):
        mat=cube[w].astype(float)                            # (512,5)
        d=well_chirality(mat)
        for k,name in enumerate(ORDER):
            A[i,k]=d[name]["gauss_auc"]; C[i,k]=d[name]["gauss_em_center"]
            Fw[i,k]=d[name]["gauss_fwhm"]; R2[i,k]=d[name].get("gauss_r2",np.nan)
    A=np.nan_to_num(A,nan=0.0); A[A<0]=0
    return A,C,Fw,R2

GROUPS={"sham":(cs,Gs["sham"]),"tbi":(ct,Gt["tbi"]),
        "resv_tbi":(ct,Gt["resv_tbi"]),"mg132_tbi":(ct,Gt["mg132_tbi"]),
        "resv_sham":(cs,Gs["resv_sham"]),"mg132_sham":(cs,Gs["mg132_sham"]),
        "swcnt_noncm":(cs,Gs["swcnt_noncm"])}
AUC={}; CEN={}; FW={}; R2={}
for g,(cube,idx) in GROUPS.items():
    AUC[g],CEN[g],FW[g],R2[g]=group_feats(cube,idx)
pr("median split-Gaussian fit R2 across wells/chiralities: %.3f"%np.nanmedian(np.vstack(list(R2.values()))))

# AUC fraction features (brightness invariant)
def frac(A): return A/np.clip(A.sum(1,keepdims=True),1e-30,None)
FF={g:frac(AUC[g]) for g in AUC}

# ---- Obj1: sham vs TBI on 12 chirality fractions ----
def loo_auc(Xa,Xb):
    X=np.vstack([Xa,Xb]); y=np.r_[np.zeros(len(Xa)),np.ones(len(Xb))]; sc=np.zeros(len(y))
    for tr,te in LeaveOneOut().split(X):
        s=StandardScaler().fit(X[tr]); l=LDA().fit(s.transform(X[tr]),y[tr]); sc[te]=l.decision_function(s.transform(X[te]))
    return roc_auc_score(y,sc),X,y
def perm(X,y,obs,n=5000,seed=0):
    rng=np.random.default_rng(seed); c=0
    for _ in range(n):
        yp=rng.permutation(y); sc=np.zeros(len(y))
        for tr,te in LeaveOneOut().split(X):
            s=StandardScaler().fit(X[tr]); l=LDA().fit(s.transform(X[tr]),yp[tr]); sc[te]=l.decision_function(s.transform(X[te]))
        if roc_auc_score(yp,sc)>=obs: c+=1
    return (c+1)/(n+1)
def cohend(a,b):
    s=np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2)); return (b.mean()-a.mean())/s
def bh(p):
    p=np.asarray(p); o=np.argsort(p); r=np.empty(len(p)); ranked=p[o]*len(p)/(np.arange(len(p))+1)
    ranked=np.minimum.accumulate(ranked[::-1])[::-1]; r[o]=np.clip(ranked,0,1); return r

auc,X,y=loo_auc(FF["sham"],FF["tbi"]); p=perm(X,y,auc,n=5000)
pr("\n[3H obj1] sham vs TBI on 12 chirality AUC-fractions: LOO-LDA AUC=%.3f perm-p=%.4f"%(auc,p))
pr("  per-chirality fraction (sham -> TBI):")
ps=[]; ds=[]
for k,name in enumerate(ORDER):
    a=FF["sham"][:,k]; b=FF["tbi"][:,k]; d=cohend(a,b); t,pp=ttest_ind(a,b); ps.append(pp); ds.append(d)
q=bh(ps)
for k,name in enumerate(ORDER):
    star="*" if q[k]<0.05 else " "
    pr("   %-8s em=%4.0fnm  sham=%.3f TBI=%.3f  d=%+.2f  p=%.1e q=%.1e %s"%(
        SHORT[name],CHIR[name]["emission_nm"],FF["sham"][:,k].mean(),FF["tbi"][:,k].mean(),ds[k],ps[k],q[k],star))
# center shifts (TBI-sham) on chiralities with decent fit
pr("  peak-center shift TBI-sham (nm):")
for k,name in enumerate(ORDER):
    cs_=np.nanmedian(CEN["sham"][:,k]); ct_=np.nanmedian(CEN["tbi"][:,k])
    pr("   %-8s %+.2f"%(SHORT[name],ct_-cs_))

# ---- Obj2: drug projection on chirality-feature axis ----
sc=StandardScaler().fit(np.vstack([FF["sham"],FF["tbi"]]))
lda=LDA().fit(sc.transform(np.vstack([FF["sham"],FF["tbi"]])),np.r_[np.zeros(len(FF["sham"])),np.ones(len(FF["tbi"]))])
flip=1.0 if lda.decision_function(sc.transform(FF["tbi"])).mean()>lda.decision_function(sc.transform(FF["sham"])).mean() else -1.0
P={g:flip*lda.decision_function(sc.transform(FF[g])) for g in FF}
pr("\n[3H obj2] drug arms on sham(0)->TBI(+) chirality-feature axis:")
for nm in ["sham","resv_sham","mg132_sham","tbi","resv_tbi","mg132_tbi","swcnt_noncm"]:
    pr("   %-12s % .3f +/- %.3f"%(nm,P[nm].mean(),P[nm].std()))
pr("   resv_tbi vs tbi d=%+.2f (neg=rescue) | mg132_tbi vs tbi d=%+.2f"%(cohend(P["tbi"],P["resv_tbi"]),cohend(P["tbi"],P["mg132_tbi"])))

# ---- figure ----
fig,ax=plt.subplots(1,3,figsize=(17,5))
# A: mean chirality fraction sham vs TBI
xk=np.arange(12); w=0.4
emc=[CHIR[n]["emission_nm"] for n in ORDER]; o=np.argsort(emc)
ms=[FF["sham"][:,k].mean() for k in o]; es=[FF["sham"][:,k].std()/np.sqrt(len(FF["sham"])) for k in o]
mt=[FF["tbi"][:,k].mean() for k in o];  et=[FF["tbi"][:,k].std()/np.sqrt(len(FF["tbi"])) for k in o]
ax[0].bar(xk-w/2,ms,w,yerr=es,color=COL["sham"],label="sham",capsize=2)
ax[0].bar(xk+w/2,mt,w,yerr=et,color=COL["tbi"],label="TBI",capsize=2)
for j,k in enumerate(o):
    if q[k]<0.05: ax[0].text(j,max(ms[j],mt[j])+0.005,"*",ha="center",fontsize=12)
ax[0].set_xticks(xk); ax[0].set_xticklabels([SHORT[ORDER[k]] for k in o],rotation=45,ha="right",fontsize=8)
ax[0].set_ylabel("chirality AUC fraction"); ax[0].set_title("3H fitted-peak chirality dist (AUC=%.3f p=%.3f)"%(auc,p),fontsize=10); ax[0].legend(fontsize=8)
# B: per-chirality effect size
ax[1].axhline(0,c="k",lw=.6); ax[1].bar(xk,[ds[k] for k in o],color=["#cb4335" if ds[k]>0 else "#2e86c1" for k in o])
for j,k in enumerate(o):
    if q[k]<0.05: ax[1].text(j,ds[k]+0.05*np.sign(ds[k] or 1),"*",ha="center",fontsize=12)
ax[1].set_xticks(xk); ax[1].set_xticklabels([SHORT[ORDER[k]] for k in o],rotation=45,ha="right",fontsize=8)
ax[1].set_ylabel("Cohen d (TBI-sham)"); ax[1].set_title("which chirality drives it",fontsize=10)
# C: drug projection
order=["sham","resv_sham","mg132_sham","tbi","resv_tbi","mg132_tbi","swcnt_noncm"]
for xi,nm in enumerate(order):
    v=P[nm]; ax[2].scatter(np.full(len(v),xi)+np.random.uniform(-.12,.12,len(v)),v,c=COL[nm],s=26,edgecolor="k",lw=.3,zorder=3)
    ax[2].plot([xi-.25,xi+.25],[v.mean(),v.mean()],c="k",lw=2,zorder=4)
ax[2].axhline(P["sham"].mean(),c=COL["sham"],ls="--",lw=1); ax[2].axhline(P["tbi"].mean(),c=COL["tbi"],ls="--",lw=1)
ax[2].set_xticks(range(len(order))); ax[2].set_xticklabels(order,rotation=40,ha="right",fontsize=8)
ax[2].set_ylabel("sham(0)->TBI(+) axis (chirality feats)"); ax[2].set_title("3H drug arms  resvD=%+.2f mgD=%+.2f"%(cohend(P["tbi"],P["resv_tbi"]),cohend(P["tbi"],P["mg132_tbi"])),fontsize=10)
fig.suptitle("exp3 3H — canonical split-Gaussian chirality-peak analysis (12 (n,m) bands)",fontsize=13)
fig.tight_layout(rect=[0,0,1,0.95])
pth=os.path.join(OUT,"exp3_3H_peakfit.png"); fig.savefig(pth,dpi=150); pr("\nsaved",pth)
