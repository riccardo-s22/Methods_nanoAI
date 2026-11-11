"""
exp3 3H — add the PAIRED chirality features (intensity RATIOS) that exp3_peakfit.py omitted.
Same canonical split-Gaussian fits (12 (n,m)); now build, per well:
  - 12 single-peak AUC fractions (as before)
  - 66 scale-invariant log10 intensity ratios log(max_A/max_B) for all chirality pairs
    (this is cie.compute_pairwise_features' ratio, the descriptor exp1's (n,m)-ratio biomarker used)
Test sham vs TBI 3H on: fractions / ratios / both, with a label-blind PCA->LOO-LDA pipeline (fast,
valid permutation null), per-ratio univariate Mann-Whitney + BH-FDR, and project the drug arms.
"""
import os, sys, warnings, importlib.util
from itertools import combinations
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
warnings.filterwarnings("ignore")
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, ".."))                    # section top level: exp3_load.py
sys.path.insert(0, os.path.join(_HERE, "..", "..", "shared"))
from repo_paths import data_path, results_dir, env_path
import numpy as np
import exp3_load as L
from scipy.stats import mannwhitneyu
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def pr(*a): print(*a); sys.stdout.flush()
OUT=env_path("EXP3_OUT", "") or results_dir("04_drug_response_prediction", "sensitivity_true_log_ratios"); os.makedirs(OUT, exist_ok=True)
CACHE_DIR=env_path("EXP3_CACHE_DIR", data_path("04_drug_response_prediction", "chirality_features"))  # shipped cache, read-only

PD=os.path.join(_HERE,"..","vendored")
spec=importlib.util.spec_from_file_location("cie",os.path.join(PD,"chirality_interpolated_extraction.py"))
cie=importlib.util.module_from_spec(spec); spec.loader.exec_module(cie)
CHIR=cie.CHIRALITY_POSITIONS; ORDER=list(CHIR.keys())
SH={k:k.replace("ch","(").replace("_",",")+")" for k in ORDER}
PAIRS=list(combinations(range(12),2))                 # 66 pairs
PNAME=[f"{SH[ORDER[a]]}/{SH[ORDER[b]]}" for a,b in PAIRS]
EXC=np.array([570,640,670,750,780],float)
COL={"sham":"#2e86c1","tbi":"#cb4335","resv_tbi":"#27ae60","mg132_tbi":"#8e44ad",
     "resv_sham":"#a9cce3","mg132_sham":"#d2b4de","swcnt_noncm":"#7f8c8d"}

cs,em=L.load_plate("SHAM","3H"); ct,_=L.load_plate("TBI","3H")
Gs=L.groups("SHAM"); Gt=L.groups("TBI"); exmin,exmax=EXC.min(),EXC.max()

def well_chir(mat):
    sp,emf,exf=cie.build_interpolated_surface(mat,em,EXC,factor=cie.INTERP_FACTOR)
    mx=np.zeros(12); au=np.zeros(12)
    for k,(name,c) in enumerate(CHIR.items()):
        tex=min(max(c["excitation_nm"],exmin),exmax)
        be,bx,m=cie.find_local_max_in_radius(sp,emf,exf,c["emission_nm"],tex,cie.SEARCH_RADIUS_NM)
        ep,ip=cie.extract_consistent_gradient_points(sp,emf,be,bx)
        d=cie.fit_gaussian_and_extract(ep,ip,be,m)
        mx[k]=max(d["gauss_max"],0); au[k]=max(np.nan_to_num(d["gauss_auc"]),0)
    return mx,au

def group(cube,idx):
    MX=np.zeros((len(idx),12)); AU=np.zeros((len(idx),12))
    for i,w in enumerate(idx): MX[i],AU[i]=well_chir(cube[w].astype(float))
    return MX,AU
GRP={"sham":(cs,Gs["sham"]),"tbi":(ct,Gt["tbi"]),"resv_tbi":(ct,Gt["resv_tbi"]),
     "mg132_tbi":(ct,Gt["mg132_tbi"]),"resv_sham":(cs,Gs["resv_sham"]),
     "mg132_sham":(cs,Gs["mg132_sham"]),"swcnt_noncm":(cs,Gs["swcnt_noncm"])}
# TRUELOG RERUN: peaks are taken from the cached canonical fits used by every reported exp3 analysis
# (refitting from raw fails in this environment: curve_fit "initial guess outside bounds").
_z=np.load(os.path.join(CACHE_DIR,"exp3_3H_peakfeats.npz"))
MX={g:_z["MX_"+g] for g in GRP}; AU={g:_z["AU_"+g] for g in GRP}

flo=float(os.environ.get("FLO","0.0"))  # TRUELOG RERUN (was 1e-6)
def fracs(A): return A/np.clip(A.sum(1,keepdims=True),1e-30,None)
def ratios(M):                                        # log10 max_A/max_B, 66 features
    R=np.zeros((len(M),len(PAIRS)))
    for j,(a,b) in enumerate(PAIRS): R[:,j]=np.log10((M[:,a]+flo)/(M[:,b]+flo))
    return R
FR={g:fracs(AU[g]) for g in AU}
RT={g:ratios(MX[g]) for g in MX}
BOTH={g:np.hstack([FR[g],RT[g]]) for g in FR}

def loo_auc(Xa,Xb,nc=15):
    X=np.vstack([Xa,Xb]); y=np.r_[np.zeros(len(Xa)),np.ones(len(Xb))]
    nc=min(nc,X.shape[1],len(y)-2)
    Z=PCA(nc,random_state=0).fit_transform(StandardScaler().fit_transform(X))
    s=np.zeros(len(y)); lda=LDA()
    for tr,te in LeaveOneOut().split(Z): lda.fit(Z[tr],y[tr]); s[te]=lda.decision_function(Z[te])
    return roc_auc_score(y,s),Z,y
def perm(Z,y,obs,n=2000,seed=0):
    rng=np.random.default_rng(seed); c=0; lda=LDA()
    for _ in range(n):
        yp=rng.permutation(y); s=np.zeros(len(y))
        for tr,te in LeaveOneOut().split(Z): lda.fit(Z[tr],yp[tr]); s[te]=lda.decision_function(Z[te])
        if roc_auc_score(yp,s)>=obs: c+=1
    return (c+1)/(n+1)
def bh(p):
    p=np.asarray(p); o=np.argsort(p); r=np.empty(len(p)); rk=p[o]*len(p)/(np.arange(len(p))+1)
    rk=np.minimum.accumulate(rk[::-1])[::-1]; r[o]=np.clip(rk,0,1); return r
def cohend(a,b):
    s=np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2)); return (b.mean()-a.mean())/s

pr("######## exp3 3H — chirality features WITH pairwise intensity ratios ########")
for nm,FEAT in [("12 AUC-fractions",FR),("66 log-ratios",RT),("both (78)",BOTH)]:
    auc,Z,y=loo_auc(FEAT["sham"],FEAT["tbi"]); p=perm(Z,y,auc,n=2000)
    pr("  [%-18s] sham vs TBI  LOO-LDA AUC=%.3f  perm-p=%.4f"%(nm,auc,p))

# univariate per-ratio sham vs TBI
pv=[]; dd=[]; aucs=[]
for j in range(len(PAIRS)):
    a=RT["sham"][:,j]; b=RT["tbi"][:,j]
    pv.append(mannwhitneyu(a,b).pvalue); dd.append(cohend(a,b))
    aucs.append(roc_auc_score(np.r_[np.zeros(len(a)),np.ones(len(b))],np.r_[a,b]))
q=bh(pv); aucs=np.array(aucs); order=np.argsort(pv)
pr("\n  top discriminative intensity ratios (sham vs TBI), BH-FDR:")
nsig=int((q<0.05).sum()); pr("   # ratios with q<0.05: %d / 66"%nsig)
for j in order[:12]:
    star="*" if q[j]<0.05 else " "
    pr("   %-16s d=%+.2f AUC=%.2f p=%.1e q=%.2f %s"%(PNAME[j],dd[j],max(aucs[j],1-aucs[j]),pv[j],q[j],star))

# drug projection on the BEST set (both), axis from sham vs TBI
def axis(FEAT):
    X=np.vstack([FEAT["sham"],FEAT["tbi"]]); y=np.r_[np.zeros(len(FEAT["sham"])),np.ones(len(FEAT["tbi"]))]
    sc=StandardScaler().fit(X); pca=PCA(min(15,X.shape[1],len(y)-2),random_state=0).fit(sc.transform(X))
    lda=LDA().fit(pca.transform(sc.transform(X)),y)
    f=1.0 if lda.decision_function(pca.transform(sc.transform(FEAT["tbi"]))).mean()>lda.decision_function(pca.transform(sc.transform(FEAT["sham"]))).mean() else -1.0
    return {g:f*lda.decision_function(pca.transform(sc.transform(FEAT[g]))) for g in FEAT}
P=axis(BOTH)
pr("\n  [drug arms on 'both'-feature sham(0)->TBI(+) axis]:")
for g in ["sham","resv_sham","mg132_sham","tbi","resv_tbi","mg132_tbi","swcnt_noncm"]:
    pr("   %-12s % .3f +/- %.3f"%(g,P[g].mean(),P[g].std()))
pr("   resv_tbi vs tbi d=%+.2f (neg=rescue) | mg132_tbi vs tbi d=%+.2f"%(cohend(P["tbi"],P["resv_tbi"]),cohend(P["tbi"],P["mg132_tbi"])))

# figure
fig,ax=plt.subplots(1,2,figsize=(14,5.5))
top=order[:14][::-1]
ax[0].barh(range(len(top)),[dd[j] for j in top],color=["#cb4335" if dd[j]>0 else "#2e86c1" for j in top])
ax[0].set_yticks(range(len(top))); ax[0].set_yticklabels([PNAME[j] for j in top],fontsize=8)
for i,j in enumerate(top):
    if q[j]<0.05: ax[0].text(dd[j]+0.02*np.sign(dd[j] or 1),i,"*",va="center",fontsize=12)
ax[0].axvline(0,c="k",lw=.6); ax[0].set_xlabel("Cohen d (TBI-sham), log10 ratio")
ax[0].set_title("top chirality intensity-RATIO differences 3H\n(%d/66 q<0.05)"%nsig,fontsize=10)
od=["sham","resv_sham","mg132_sham","tbi","resv_tbi","mg132_tbi","swcnt_noncm"]
for xi,g in enumerate(od):
    v=P[g]; ax[1].scatter(np.full(len(v),xi)+np.random.uniform(-.12,.12,len(v)),v,c=COL[g],s=26,edgecolor="k",lw=.3,zorder=3)
    ax[1].plot([xi-.25,xi+.25],[v.mean(),v.mean()],c="k",lw=2,zorder=4)
ax[1].axhline(P["sham"].mean(),c=COL["sham"],ls="--",lw=1); ax[1].axhline(P["tbi"].mean(),c=COL["tbi"],ls="--",lw=1)
ax[1].set_xticks(range(len(od))); ax[1].set_xticklabels(od,rotation=40,ha="right",fontsize=8)
ax[1].set_ylabel("sham(0)->TBI(+) axis (fractions+ratios)"); ax[1].set_title("3H drug arms on paired-feature axis",fontsize=10)
fig.suptitle("exp3 3H — fitted-peak chirality features WITH pairwise intensity ratios",fontsize=12)
fig.tight_layout(rect=[0,0,1,0.95])
pth=os.path.join(OUT,"exp3_3H_peakfit_ratios.png"); fig.savefig(pth,dpi=150); pr("\nsaved",pth)
