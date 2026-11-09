"""
exp3 3H — four drug-effect analyses exploiting the 2x3 factorial (condition{sham,TBI} x
drug{none,resveratrol,MG132}; each drug on BOTH plates). Features = cached canonical fitted-peak
12 fractions + 66 log-ratios.

(1) DiD / interaction: rescue = (resv_TBI - TBI) - (resv_sham - sham); cancels additive plate & drug
    main effects. On the sham->TBI axis biomarker and the interpretable (6,5)/(7,6) ratio; also an
    OLS with well-column covariate.
(2) Within-plate PERMANOVA on the TBI plate only (TBI vs resv-TBI vs MG132-TBI) -> batch-free.
(3) Equivalence (TOST): is resv-TBI equivalent to sham? + % rescue with bootstrap CI.
(4) Per-chirality rescue map: DiD interaction per feature, BH-FDR -> which drivers each drug normalizes.
"""
import os, sys, warnings
from itertools import combinations
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
warnings.filterwarnings("ignore")
import numpy as np
import exp3_load as L
from scipy import stats
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.preprocessing import StandardScaler
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def pr(*a): print(*a); sys.stdout.flush()
OUT="/mnt/c/Users/riccardo-s/Documents/CNT/TBI/exp3"; TP=os.environ.get("EXP3_TP","3H")
ORDER=["ch8_3","ch6_5","ch7_5","ch10_2","ch9_4","ch8_4","ch7_6","ch8_6","ch8_7","ch9_5","ch10_3","ch10_5"]
SHn={k:k.replace("ch","(").replace("_",",")+")" for k in ORDER}
PAIRS=list(combinations(range(12),2)); PNAME=[f"{SHn[ORDER[a]]}/{SHn[ORDER[b]]}" for a,b in PAIRS]
FEAT=[f"frac{SHn[ORDER[k]]}" for k in range(12)]+[f"r:{n}" for n in PNAME]
flo=1e-6
R65_76=12+PAIRS.index((ORDER.index("ch6_5"),ORDER.index("ch7_6")))   # feature idx of (6,5)/(7,6)

z=np.load("exp3_%s_peakfeats.npz"%TP)
def feats(g):
    A=z[f"AU_{g}"]; M=z[f"MX_{g}"]; fr=A/np.clip(A.sum(1,keepdims=True),1e-30,None)
    Rr=np.column_stack([np.log10((M[:,a]+flo)/(M[:,b]+flo)) for a,b in PAIRS]); return np.hstack([fr,Rr])
GR=["sham","tbi","resv_tbi","mg132_tbi","resv_sham","mg132_sham"]
X={g:feats(g) for g in GR}
# well columns (for covariate); from plate layout
Gs=L.groups("SHAM"); Gt=L.groups("TBI")
IDX={"sham":Gs["sham"],"resv_sham":Gs["resv_sham"],"mg132_sham":Gs["mg132_sham"],
     "tbi":Gt["tbi"],"resv_tbi":Gt["resv_tbi"],"mg132_tbi":Gt["mg132_tbi"]}
COLn={g:(IDX[g]%12+1).astype(float) for g in GR}

# ----- biomarkers -----
# axis: sham vs TBI LDA on 66 ratios, oriented TBI+
def ratio_only(M): return M[:,12:]
sc=StandardScaler().fit(np.vstack([ratio_only(X["sham"]),ratio_only(X["tbi"])]))
lda=LDA().fit(sc.transform(np.vstack([ratio_only(X["sham"]),ratio_only(X["tbi"])])),
              np.r_[np.zeros(len(X["sham"])),np.ones(len(X["tbi"]))])
flip=1.0 if lda.decision_function(sc.transform(ratio_only(X["tbi"]))).mean()>lda.decision_function(sc.transform(ratio_only(X["sham"]))).mean() else -1.0
AX={g:flip*lda.decision_function(sc.transform(ratio_only(X[g]))) for g in GR}    # sham<low TBI>high
R6576={g:X[g][:,R65_76] for g in GR}                                              # interpretable ratio

def did_contrast(bm, drug_tbi, drug_sham):
    """ (drug_TBI - TBI) - (drug_sham - sham); SE from 4 group variances. <0 = rescue (axis TBI+)."""
    a,b=bm["tbi"],bm[drug_tbi]; c,d=bm["sham"],bm[drug_sham]
    did=(b.mean()-a.mean())-(d.mean()-c.mean())
    se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b)+c.var(ddof=1)/len(c)+d.var(ddof=1)/len(d))
    t=did/se; df=len(a)+len(b)+len(c)+len(d)-4; p=2*stats.t.sf(abs(t),df)
    return did,se,t,p

pr("############ exp3 %s drug effects (2x3 factorial) ############"%TP)
pr("\n=== (1) DIFFERENCE-IN-DIFFERENCES / interaction (rescue) ===")
pr(" biomarker: sham->TBI ratio-axis (TBI+).  raw group means:")
for g in GR: pr("   %-12s % .3f +/- %.3f"%(g,AX[g].mean(),AX[g].std()))
for bmname,bm in [("axis",AX),("(6,5)/(7,6)",R6576)]:
    pr(" [%s] DiD interaction (neg = rescue toward sham):"%bmname)
    for label,dt,ds in [("resveratrol","resv_tbi","resv_sham"),("MG132","mg132_tbi","mg132_sham")]:
        did,se,t,p=did_contrast(bm,dt,ds)
        simple=bm[dt].mean()-bm["tbi"].mean()   # naive within-TBI-plate effect (no DiD)
        pr("   %-12s DiD=%+.3f (SE %.3f) t=%+.2f p=%.3f | naive(%s-TBI)=%+.3f"%(label,did,se,t,p,label[:4],simple))
# OLS with column covariate (axis biomarker): y ~ cond + drug2 + drug3 + cond:drug2 + cond:drug3 + col
pr(" [axis] OLS interaction with well-column covariate:")
yv=[]; cond=[]; dr=[]; col=[]
for g in GR:
    yv.append(AX[g]); cond.append(np.full(len(AX[g]),0.0 if "sham" in g and g!="mg132_sham" else 0.0))
yv=np.concatenate([AX[g] for g in GR])
condv=np.concatenate([np.full(len(AX[g]),1.0 if g in ("tbi","resv_tbi","mg132_tbi") else 0.0) for g in GR])
d_resv=np.concatenate([np.full(len(AX[g]),1.0 if "resv" in g else 0.0) for g in GR])
d_mg=np.concatenate([np.full(len(AX[g]),1.0 if "mg132" in g else 0.0) for g in GR])
colv=np.concatenate([COLn[g] for g in GR])
Dm=np.column_stack([np.ones_like(yv),condv,d_resv,d_mg,condv*d_resv,condv*d_mg,colv])
beta,_,_,_=np.linalg.lstsq(Dm,yv,rcond=None); resid=yv-Dm@beta
dof=len(yv)-Dm.shape[1]; s2=resid@resid/dof; cov=s2*np.linalg.inv(Dm.T@Dm); seb=np.sqrt(np.diag(cov))
names=["intercept","condition(TBI)","drug(resv)","drug(MG132)","TBI:resv","TBI:MG132","well_col"]
for i,nm in enumerate(names):
    tt=beta[i]/seb[i]; pp=2*stats.t.sf(abs(tt),dof)
    star="*" if pp<0.05 else " "
    pr("   %-16s coef=%+.3f SE=%.3f t=%+.2f p=%.3f %s"%(nm,beta[i],seb[i],tt,pp,star))

# ----- (2) within-plate PERMANOVA (TBI plate) -----
def permanova(groups, B=5000, seed=0):
    Xs=np.vstack([groups[g] for g in groups]); lab=np.concatenate([[g]*len(groups[g]) for g in groups])
    Xz=StandardScaler().fit_transform(Xs); D=np.linalg.norm(Xz[:,None,:]-Xz[None,:,:],axis=2)**2
    N=len(lab); SST=D.sum()/(2*N)
    def ssw(l):
        s=0.0
        for g in set(l):
            m=np.where(l==g)[0]; s+=D[np.ix_(m,m)].sum()/(2*len(m))
        return s
    a=len(set(lab)); SSW=ssw(lab); F=((SST-SSW)/(a-1))/(SSW/(N-a))
    rng=np.random.default_rng(seed); c=0
    for _ in range(B):
        lp=rng.permutation(lab); SWp=ssw(lp); Fp=((SST-SWp)/(a-1))/(SWp/(N-a))
        if Fp>=F: c+=1
    return F,(c+1)/(B+1)
pr("\n=== (2) WITHIN-PLATE PERMANOVA (batch-free; 78 chirality feats, standardized) ===")
F,p=permanova({g:X[g] for g in ["tbi","resv_tbi","mg132_tbi"]})
pr("   TBI plate 3-group (TBI/resv/MG132): pseudo-F=%.2f perm-p=%.4f"%(F,p))
for pair in [("tbi","resv_tbi"),("tbi","mg132_tbi"),("resv_tbi","mg132_tbi")]:
    F,p=permanova({g:X[g] for g in pair}); pr("     %-22s F=%.2f p=%.4f"%(" vs ".join(pair),F,p))
F,p=permanova({g:X[g] for g in ["sham","resv_sham","mg132_sham"]})
pr("   SHAM plate 3-group (ref): pseudo-F=%.2f perm-p=%.4f"%(F,p))

# ----- (3) equivalence / TOST + % rescue -----
def tost(a,b,margin):
    # H0: |mean diff| >= margin ; equivalence if both one-sided tests reject
    diff=b.mean()-a.mean(); se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b)); df=len(a)+len(b)-2
    t1=(diff+margin)/se; p1=stats.t.sf(t1,df)      # H0: diff<=-margin
    t2=(diff-margin)/se; p2=stats.t.cdf(t2,df)     # H0: diff>=+margin
    return diff,max(p1,p2)
pr("\n=== (3) EQUIVALENCE (TOST): is resv-TBI equivalent to sham? ===")
gap=AX["tbi"].mean()-AX["sham"].mean(); sd_sham=AX["sham"].std(ddof=1)
for margfac,desc in [(0.5,"0.5*sham_sd"),(0.5*gap/sd_sham,"0.5*(TBI-sham)gap")]:
    marg=margfac*sd_sham
    diff,pe=tost(AX["sham"],AX["resv_tbi"],marg)
    pr("   resv_TBI vs sham: diff=%+.3f, margin=%.3f (%s) -> TOST p=%.3f %s"%(diff,marg,desc,pe,"EQUIV" if pe<0.05 else "not equiv"))
# % rescue with bootstrap CI
def pct_rescue(seed=0):
    rng=np.random.default_rng(seed); boots=[]
    for _ in range(5000):
        s=AX["sham"][rng.integers(0,len(AX["sham"]),len(AX["sham"]))].mean()
        t=AX["tbi"][rng.integers(0,len(AX["tbi"]),len(AX["tbi"]))].mean()
        r=AX["resv_tbi"][rng.integers(0,len(AX["resv_tbi"]),len(AX["resv_tbi"]))].mean()
        if abs(t-s)>1e-9: boots.append(100*(t-r)/(t-s))
    return np.median(boots),np.percentile(boots,2.5),np.percentile(boots,97.5)
m,lo,hi=pct_rescue()
pr("   resveratrol %% rescue toward sham = %.0f%% (95%% CI %.0f..%.0f%%)  [100%%=full sham, 0%%=still TBI]"%(m,lo,hi))

# ----- (4) per-chirality rescue map -----
pr("\n=== (4) PER-CHIRALITY rescue map: DiD interaction per feature (resveratrol), BH-FDR ===")
def bh(p):
    p=np.asarray(p); o=np.argsort(p); r=np.empty(len(p)); rk=p[o]*len(p)/(np.arange(len(p))+1)
    rk=np.minimum.accumulate(rk[::-1])[::-1]; r[o]=np.clip(rk,0,1); return r
did=np.zeros(78); pp=np.zeros(78)
for j in range(78):
    bm={g:X[g][:,j] for g in GR}; d,se,t,p=did_contrast(bm,"resv_tbi","resv_sham"); did[j]=d; pp[j]=p
q=bh(pp); od=np.argsort(pp)
# orient sign by TBI direction so negative DiD = normalizing
dvec=np.array([X["tbi"][:,j].mean()-X["sham"][:,j].mean() for j in range(78)])
norm_score=did*np.sign(dvec)   # <0 => moves feature back toward sham (rescue)
pr("   top features by DiD p (resv); norm<0 = resveratrol normalizes toward sham:")
for j in od[:10]:
    pr("     %-18s DiD=%+.3f norm=%+.3f p=%.3f q=%.2f"%(FEAT[j],did[j],norm_score[j],pp[j],q[j]))
pr("   #features q<0.05: %d/78 ; net normalization (mean norm_score over TBI-driver feats |d|>0.3): %.3f"
   %(int((q<0.05).sum()),norm_score[np.abs(dvec)>0.3*np.std(dvec)].mean()))

# ----- figure -----
fig,ax=plt.subplots(1,3,figsize=(17,5.2))
order6=["sham","resv_sham","mg132_sham","tbi","resv_tbi","mg132_tbi"]
colc={"sham":"#2e86c1","resv_sham":"#a9cce3","mg132_sham":"#d2b4de","tbi":"#cb4335","resv_tbi":"#27ae60","mg132_tbi":"#8e44ad"}
for xi,g in enumerate(order6):
    v=AX[g]; ax[0].scatter(np.full(len(v),xi)+np.random.uniform(-.12,.12,len(v)),v,c=colc[g],s=24,edgecolor="k",lw=.3,zorder=3)
    ax[0].plot([xi-.25,xi+.25],[v.mean()]*2,c="k",lw=2,zorder=4)
ax[0].axhline(AX["sham"].mean(),c="#2e86c1",ls="--",lw=1); ax[0].axhline(AX["tbi"].mean(),c="#cb4335",ls="--",lw=1)
ax[0].set_xticks(range(6)); ax[0].set_xticklabels(order6,rotation=40,ha="right",fontsize=8)
ax[0].set_ylabel("sham->TBI axis"); ax[0].set_title("(1) factorial on axis biomarker",fontsize=10)
# DiD bars
labs=["resv","MG132"]; dv=[did_contrast(AX,"resv_tbi","resv_sham")[0],did_contrast(AX,"mg132_tbi","mg132_sham")[0]]
ev=[did_contrast(AX,"resv_tbi","resv_sham")[1],did_contrast(AX,"mg132_tbi","mg132_sham")[1]]
ax[1].bar(labs,dv,yerr=ev,color=["#27ae60","#8e44ad"],capsize=4); ax[1].axhline(0,c="k",lw=.6)
ax[1].set_ylabel("DiD interaction (neg=rescue)"); ax[1].set_title("(1) DiD: drug x condition",fontsize=10)
# per-chirality rescue map (top normalized)
top=od[:14]; top=top[np.argsort([norm_score[j] for j in top])]
ax[2].barh(range(len(top)),[norm_score[j] for j in top],color=["#27ae60" if norm_score[j]<0 else "#cb4335" for j in top])
ax[2].set_yticks(range(len(top))); ax[2].set_yticklabels([FEAT[j] for j in top],fontsize=7); ax[2].axvline(0,c="k",lw=.6)
ax[2].set_xlabel("resv DiD x TBI-sign (neg=normalized)"); ax[2].set_title("(4) per-chirality rescue map",fontsize=10)
fig.suptitle("exp3 %s — drug-effect analyses (DiD interaction / PERMANOVA / TOST / per-chirality)"%TP,fontsize=12)
fig.tight_layout(rect=[0,0,1,0.95]); p=os.path.join(OUT,"exp3_%s_drug_effects.png"%TP); fig.savefig(p,dpi=150); pr("\nsaved",p)
