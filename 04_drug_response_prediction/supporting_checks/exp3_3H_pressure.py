"""
exp3 3H — PRESSURE-TEST the headline sham-vs-TBI result (shape AUC~0.90, Gaussian-ratio AUC~0.94)
against spatial gradient and outlier wells.

Design fact: sham = SHAM-plate cols1-4 (all 8 rows); tbi = TBI-plate cols5-8 (all 8 rows).
 -> ROW gradient CANNOT fake sham-vs-TBI (identical row span in both groups).
 -> Only a COLUMN gradient could contribute (cols1-4 vs cols5-8), and it is confounded with
    plate+treatment, so it cannot be fully separated — but we can QUANTIFY it and bound its effect.

Battery:
 [A] per-plate intensity gradient (row vs col slopes, R2) — is there a COLUMN gradient at all?
 [B] outlier wells per plate (L.outlier_keep) -> re-run AUC (shape) with outliers removed
 [C] position residualization: remove within-group [row,col] trend -> re-run AUC
 [D] alignment: cosine(sham->TBI shape-diff, per-plate spatial-gradient directions) &
     fraction of the sham->TBI difference explained by a same-treatment COLUMN-gradient estimate
 [E] Gaussian-ratio AUC re-run with outliers removed (headline 0.94)
"""
import os, sys, warnings, importlib.util
from itertools import combinations
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
warnings.filterwarnings("ignore")
import numpy as np
if not hasattr(np,"trapz"): np.trapz=np.trapezoid
import exp3_load as L
import exp3_7d_load as L7            # reuse outlier_keep
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score
def pr(*a): print(*a); sys.stdout.flush()
OUT="/mnt/c/Users/riccardo-s/Documents/CNT/TBI/exp3"

cs,em=L.load_plate("SHAM","3H"); Gs=L.groups("SHAM")
ct,_ =L.load_plate("TBI","3H");  Gt=L.groups("TBI")

def shape_feats(cube, idx):
    X=cube[idx].astype(float); X=X/X.sum(1,keepdims=True); return X.reshape(len(idx),-1)
def loo_auc(Xa,Xb,ncomp=20):
    X=np.vstack([Xa,Xb]); y=np.r_[np.zeros(len(Xa)),np.ones(len(Xb))]
    nc=min(ncomp,X.shape[1],len(y)-2)
    Z=PCA(nc,random_state=0).fit_transform(StandardScaler().fit_transform(X))
    s=np.zeros(len(y)); lda=LDA()
    for tr,te in LeaveOneOut().split(Z): lda.fit(Z[tr],y[tr]); s[te]=lda.decision_function(Z[te])
    return roc_auc_score(y,s)

def plate_grad(cube, name):
    N=cube.shape[0]; row=np.array([i//12 for i in range(N)],float); col=np.array([i%12 for i in range(N)],float)
    tot=cube.sum(axis=(1,2)); pos=tot[tot>0].min()
    ltot=np.log(np.clip(tot,pos,None))
    def fit(cols):
        Xd=np.column_stack([np.ones(N)]+cols); b,*_=np.linalg.lstsq(Xd,ltot,rcond=None)
        r2=1-((ltot-Xd@b)**2).sum()/((ltot-ltot.mean())**2).sum(); return b,r2
    b,r2=fit([row,col]); _,r2r=fit([row]); _,r2c=fit([col])
    pr("  [%s] %+.1f%%/row  %+.1f%%/col  R2(row+col)=%.3f  R2(row only)=%.3f  R2(col only)=%.3f"
       %(name,100*(np.exp(b[1])-1),100*(np.exp(b[2])-1),r2,r2r,r2c))
    return b

pr("################ exp3 3H — PRESSURE TEST (gradient + outliers) ################")

# ---- reproduce headline ----
Xs=shape_feats(cs,Gs["sham"]); Xt=shape_feats(ct,Gt["tbi"])
auc0=loo_auc(Xs,Xt)
pr("\n[headline] sham(n=%d) vs TBI(n=%d) shape LOO-LDA AUC = %.3f  (published ~0.903)"%(len(Xs),len(Xt),auc0))

# ---- [A] gradient per plate ----
pr("\n[A] per-plate intensity gradient (log-intensity trend):")
plate_grad(cs,"SHAM"); plate_grad(ct,"TBI")

# ---- [B] outliers ----
keep_s=L7.outlier_keep(cs); keep_t=L7.outlier_keep(ct)
sham_c=np.array([i for i in Gs["sham"] if keep_s[i]]); tbi_c=np.array([i for i in Gt["tbi"] if keep_t[i]])
pr("\n[B] outlier wells removed: SHAM %d/%d wells kept (sham grp %d->%d) | TBI %d/%d (tbi grp %d->%d)"
   %(int(keep_s.sum()),len(keep_s),len(Gs["sham"]),len(sham_c),int(keep_t.sum()),len(keep_t),len(Gt["tbi"]),len(tbi_c)))
def _wl(idx): return " ".join("%s%d"%("ABCDEFGH"[i//12],i%12+1) for i in idx)
pr("    removed from sham group: %s"%(_wl([i for i in Gs['sham'] if not keep_s[i]]) or "(none)"))
pr("    removed from tbi  group: %s"%(_wl([i for i in Gt['tbi'] if not keep_t[i]]) or "(none)"))
auc_out=loo_auc(shape_feats(cs,sham_c),shape_feats(ct,tbi_c))
pr("    sham vs TBI shape AUC (outliers removed) = %.3f"%auc_out)

# ---- [C] position residualization (remove within-group row+col trend) ----
def detrend(cube, idx):
    X=shape_feats(cube,idx); r=np.array([i//12 for i in idx],float); c=np.array([i%12 for i in idx],float)
    Xd=np.column_stack([np.ones(len(idx)),r,c]); B=np.linalg.pinv(Xd)@X
    return X - Xd@B + X.mean(0)     # remove linear position trend, keep group mean
auc_dt=loo_auc(detrend(cs,Gs["sham"]),detrend(ct,Gt["tbi"]))
auc_dt_out=loo_auc(detrend(cs,sham_c),detrend(ct,tbi_c))
pr("\n[C] sham vs TBI shape AUC (within-group row+col detrended)      = %.3f"%auc_dt)
pr("    sham vs TBI shape AUC (detrended AND outliers removed)       = %.3f"%auc_dt_out)

# ---- [D] alignment of the sham->TBI axis with spatial-gradient directions ----
def cos(a,b): return float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30))
vdiff=Xt.mean(0)-Xs.mean(0)                       # sham->TBI shape difference (2560,)
# same-treatment gradient direction in feature space per plate:
# demean each treatment group, regress residual features on [row,col] -> col/row coefficient vectors
def grad_dirs(cube, groupsdict, names):
    N=cube.shape[0]; feats=[]; rr=[]; ccv=[]
    for nm in names:
        idx=groupsdict[nm]; X=shape_feats(cube,idx); X=X-X.mean(0)   # demean treatment
        feats.append(X); rr.append(np.array([i//12 for i in idx],float)-np.mean([i//12 for i in idx]))
        ccv.append(np.array([i%12 for i in idx],float)-np.mean([i%12 for i in idx]))
    X=np.vstack(feats); r=np.concatenate(rr); c=np.concatenate(ccv)
    Xd=np.column_stack([r,c]); B=np.linalg.pinv(Xd)@X     # (2, 2560): [row_dir; col_dir]
    return B[0],B[1]
grow_s,gcol_s=grad_dirs(cs,Gs,["sham","resv_sham","mg132_sham"])
grow_t,gcol_t=grad_dirs(ct,Gt,["mg132_tbi","tbi","resv_tbi"])
pr("\n[D] alignment of sham->TBI axis with SAME-TREATMENT spatial gradients (cosine):")
pr("    cos(diff, SHAM col-grad)=%+.2f  cos(diff, SHAM row-grad)=%+.2f"%(cos(vdiff,gcol_s),cos(vdiff,grow_s)))
pr("    cos(diff, TBI  col-grad)=%+.2f  cos(diff, TBI  row-grad)=%+.2f"%(cos(vdiff,gcol_t),cos(vdiff,grow_t)))
# how much shape offset would the COLUMN gradient predict between col-pos 1.5 (sham) and 5.5 (tbi)?
dcol=(5.5-1.5)  # 0-based mean col: sham cols0-3=1.5, tbi cols4-7=5.5
pred_s=gcol_s*dcol; pred_t=gcol_t*dcol
pr("    |col-grad predicted offset (SHAM slope)| / |observed sham->TBI diff| = %.2f  (cos=%+.2f)"
   %(np.linalg.norm(pred_s)/np.linalg.norm(vdiff),cos(pred_s,vdiff)))
pr("    |col-grad predicted offset (TBI  slope)| / |observed sham->TBI diff| = %.2f  (cos=%+.2f)"
   %(np.linalg.norm(pred_t)/np.linalg.norm(vdiff),cos(pred_t,vdiff)))

# ---- [E] Gaussian-ratio AUC with outliers removed ----
PD="/mnt/c/Users/riccardo-s/Documents/CNT/targetALS/exp5/physical_descriptors"
spec=importlib.util.spec_from_file_location("cie",os.path.join(PD,"chirality_interpolated_extraction.py"))
cie=importlib.util.module_from_spec(spec); spec.loader.exec_module(cie)
CHIR=cie.CHIRALITY_POSITIONS; EXCa=np.array(L.EXC,float); PAIRS=list(combinations(range(12),2))
def well_chir(mat):
    sp,emf,exf=cie.build_interpolated_surface(mat,em,EXCa,factor=cie.INTERP_FACTOR)
    mx=np.zeros(12)
    for k,(nm,cc) in enumerate(CHIR.items()):
        tex=min(max(cc["excitation_nm"],EXCa.min()),EXCa.max())
        be,bx,m=cie.find_local_max_in_radius(sp,emf,exf,cc["emission_nm"],tex,cie.SEARCH_RADIUS_NM)
        ep,ip=cie.extract_consistent_gradient_points(sp,emf,be,bx)
        d=cie.fit_gaussian_and_extract(ep,ip,be,m); mx[k]=max(d["gauss_max"],0)
    return mx
def ratios(cube,idx):
    M=np.array([well_chir(cube[w].astype(float)) for w in idx]); flo=1e-6
    R=np.zeros((len(idx),len(PAIRS)))
    for j,(a,b) in enumerate(PAIRS): R[:,j]=np.log10((M[:,a]+flo)/(M[:,b]+flo))
    return R
def loo_auc_pca(Xa,Xb,nc=15):
    X=np.vstack([Xa,Xb]); y=np.r_[np.zeros(len(Xa)),np.ones(len(Xb))]; nc=min(nc,X.shape[1],len(y)-2)
    Z=PCA(nc,random_state=0).fit_transform(StandardScaler().fit_transform(X))
    s=np.zeros(len(y)); lda=LDA()
    for tr,te in LeaveOneOut().split(Z): lda.fit(Z[tr],y[tr]); s[te]=lda.decision_function(Z[te])
    return roc_auc_score(y,s)
Rs=ratios(cs,Gs["sham"]); Rt=ratios(ct,Gt["tbi"])
auc_r=loo_auc_pca(Rs,Rt)
Rs_c=ratios(cs,sham_c); Rt_c=ratios(ct,tbi_c)
auc_r_out=loo_auc_pca(Rs_c,Rt_c)
pr("\n[E] Gaussian-ratio (66 log-ratio) sham vs TBI AUC:  full=%.3f (published ~0.938) | outliers-removed=%.3f"
   %(auc_r,auc_r_out))

pr("\n================ SUMMARY ================")
pr(" shape AUC: headline=%.3f | outliers-removed=%.3f | detrended=%.3f | detrended+clean=%.3f"%(auc0,auc_out,auc_dt,auc_dt_out))
pr(" Gaussian-ratio AUC: full=%.3f | outliers-removed=%.3f"%(auc_r,auc_r_out))
