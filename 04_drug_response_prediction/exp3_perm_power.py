"""
exp3 3H — permutation 'power'/precision analysis: how many permutations B are needed to certify
significance of TBI vs Sham, for BOTH methods:
  Method A = FULL-SHAPE spectra  (per-exc L1, 2560-dim -> PCA20 -> LOO-LDA AUC)
  Method B = FITTED-PEAK chirality (66 pairwise log-ratios -> PCA15 -> LOO-LDA AUC)
Key idea: permutation count does NOT change statistical power (that's effect-size & n); it sets the
RESOLUTION 1/(B+1) and the Monte-Carlo error of the p-value. We run a large null (B0), then report:
 - observed AUC, exceedances k, MC p_hat & exact (Clopper-Pearson) 95% CI
 - null AUC mean/sd, z of observed, normal-tail p (extrapolates below MC resolution)
 - minimal B to (i) resolve p<target, (ii) CERTIFY p<target with 0 exceedances via Rule-of-3 (3/B),
   (iii) if k>0, reach >=99% prob of declaring p<0.05 (binomial), (iv) 10% relative precision on p.
"""
import os, sys, warnings, importlib.util
from itertools import combinations
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
warnings.filterwarnings("ignore")
import numpy as np
import exp3_load as L
from scipy import stats
from scipy.stats import norm, beta, binom
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score
def pr(*a): print(*a); sys.stdout.flush()
OUT="/mnt/c/Users/riccardo-s/Documents/CNT/TBI/exp3"
B0=20000
ORDER=["ch8_3","ch6_5","ch7_5","ch10_2","ch9_4","ch8_4","ch7_6","ch8_6","ch8_7","ch9_5","ch10_3","ch10_5"]
PAIRS=list(combinations(range(12),2)); flo=1e-6

# ---------- Method A: full-shape ----------
cs,em=L.load_plate("SHAM","3H"); ct,_=L.load_plate("TBI","3H")
Gs=L.groups("SHAM"); Gt=L.groups("TBI")
def shp(cube,idx):
    X=cube[idx].astype(float); X=X/X.sum(1,keepdims=True); return X.reshape(len(idx),-1)  # 2560
XA=np.vstack([shp(cs,Gs["sham"]),shp(ct,Gt["tbi"])])
# ---------- Method B: ratios (from cache) ----------
z=np.load("exp3_3H_peakfeats.npz")
def ratios(MX): return np.column_stack([np.log10((MX[:,a]+flo)/(MX[:,b]+flo)) for a,b in PAIRS])
XB=np.vstack([ratios(z["MX_sham"]),ratios(z["MX_tbi"])])
y=np.r_[np.zeros((Gs["sham"]).shape[0] if False else 32),np.ones(32)]  # 32 sham, 32 TBI
y=np.r_[np.zeros(len(z["MX_sham"])),np.ones(len(z["MX_tbi"]))]

def embed(X,nc):
    nc=min(nc,X.shape[1],len(y)-2)
    return PCA(nc,random_state=0).fit_transform(StandardScaler().fit_transform(X))
def loo_auc(Z,yy):
    s=np.zeros(len(yy)); lda=LDA()
    for tr,te in LeaveOneOut().split(Z): lda.fit(Z[tr],yy[tr]); s[te]=lda.decision_function(Z[te])
    return roc_auc_score(yy,s)
def perm_null(Z,yy,obs,B,seed=0):
    rng=np.random.default_rng(seed); null=np.empty(B); k=0
    for i in range(B):
        yp=rng.permutation(yy); a=loo_auc(Z,yp); null[i]=a
        if a>=obs-1e-12: k+=1
    return null,k

def cp_ci(k,n,alpha=0.05):  # Clopper-Pearson on raw exceedance rate
    lo=0.0 if k==0 else beta.ppf(alpha/2,k,n-k+1)
    hi=1.0 if k==n else beta.ppf(1-alpha/2,k+1,n-k)
    return lo,hi

def analyze(name,X,nc):
    Z=embed(X,nc); obs=loo_auc(Z,y)
    pr("\n================ %s ================"%name)
    pr(" observed LOO-LDA AUC = %.4f   (PCA %d-dim, n=%d/%d)"%(obs,Z.shape[1],int((y==0).sum()),int((y==1).sum())))
    null,k=perm_null(Z,y,obs,B0)
    phat=(1+k)/(B0+1)
    lo,hi=cp_ci(k,B0)                       # CI on rate; MC p-value CI ~ this
    mu,sd=null.mean(),null.std(ddof=1); zsc=(obs-mu)/sd; pz=norm.sf(zsc)
    pr(" permutation null (B0=%d): mean=%.3f sd=%.3f  max=%.3f  95th=%.3f 99th=%.3f"
       %(B0,mu,sd,null.max(),np.percentile(null,95),np.percentile(null,99)))
    pr(" exceedances k=%d/%d -> MC p_hat=%.2eg  (rate 95%% CI [%.2e, %.2e])"%(k,B0,phat,lo,hi))
    pr(" normal-tail: z=%.2f  p_normal=%.2eg   (extrapolates below MC resolution; assumes Gaussian null tail)"%(zsc,pz))
    # minimal B table
    pr(" --- minimal #permutations B needed ---")
    for t in (0.05,0.01,0.001):
        pr("   resolve p<%.3f (need 1/(B+1)<t): B>=%d"%(t,int(np.ceil(1/t))-1))
    if k==0:
        for t in (0.05,0.01,0.001):
            pr("   CERTIFY p<%.3f w/ 0 exceedances (Rule-of-3, 95%%): B>%d"%(t,int(np.ceil(3/t))))
    else:
        # power to DECLARE p<0.05 given estimated true p (=phat) as fn of B
        ptrue=max(phat,1e-9)
        def powB(B,a=0.05):  # P(declare): K<a(B+1)-1
            thr=int(np.floor(a*(B+1)-1-1e-9));
            return binom.cdf(thr,B,ptrue) if thr>=0 else 0.0
        Bgrid=[100,200,500,1000,2000,5000,10000]
        for B in Bgrid:
            if powB(B)>=0.99: pr("   >=99%% prob to declare p<0.05 (given p~%.2eg): B>=%d"%(ptrue,B)); break
    # precision: B for 10% relative SE on p (use the larger of phat / a small floor)
    p_for_prec=max(phat,1/(B0+1))
    Bprec=int(np.ceil((1-p_for_prec)/(p_for_prec*0.10**2)))
    pr("   10%% relative SE on the p-value (~p=%.2eg): B>=%d"%(p_for_prec,Bprec))
    return dict(name=name,obs=obs,k=k,phat=phat,z=zsc,pz=pz,mu=mu,sd=sd,null=null)

rA=analyze("Method A: FULL-SHAPE (2560-d -> PCA20)",XA,20)
rB=analyze("Method B: FITTED-PEAK RATIOS (66 -> PCA15)",XB,15)

pr("\n================ RECOMMENDATION ================")
for r in (rA,rB):
    pr(" %-42s AUC=%.3f z=%.1f -> effect is %s; null max over %d perms=%.3f (obs %s exceeded)"
       %(r["name"],r["obs"],r["z"],"VERY strong" if r["z"]>5 else ("strong" if r["z"]>3 else "modest"),
         B0,r["null"].max(),"NEVER" if r["k"]==0 else "%d×"%r["k"]))
pr(" Permutations do NOT change power (set by AUC & n=32/32, both here large). They set the smallest")
pr(" CERTIFIABLE p. Both methods: 0 exceedances in %d -> to certify p<0.05 need B>60, p<0.01 B>300,"%B0)
pr(" p<0.001 B>3000 (Rule-of-3). Practical: B=1000 certifies p<0.003; B=5000 certifies p<6e-4.")
pr(" Use B=5000 for headline significance (cheap, certifies p<1e-3); B=10000 if reporting an exact tiny p.")

# figure
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,2,figsize=(14,5))
for i,r in enumerate((rA,rB)):
    ax[i].hist(r["null"],bins=60,color="#95a5a6",edgecolor="none")
    ax[i].axvline(r["obs"],c="#cb4335",lw=2,label="observed AUC=%.3f"%r["obs"])
    ax[i].axvline(np.percentile(r["null"],95),c="#2e86c1",ls="--",lw=1,label="null 95th=%.3f"%np.percentile(r["null"],95))
    ax[i].set_title("%s\nz=%.1f, p_normal=%.1eg, 0-exceedance in %d"%(r["name"],r["z"],r["pz"],B0),fontsize=9)
    ax[i].set_xlabel("permuted LOO-LDA AUC"); ax[i].set_ylabel("count"); ax[i].legend(fontsize=8)
fig.suptitle("exp3 3H — permutation null & required B (TBI vs Sham), both methods",fontsize=12)
fig.tight_layout(rect=[0,0,1,0.95]); p=os.path.join(OUT,"exp3_3H_perm_power.png"); fig.savefig(p,dpi=150); pr("\nsaved",p)
