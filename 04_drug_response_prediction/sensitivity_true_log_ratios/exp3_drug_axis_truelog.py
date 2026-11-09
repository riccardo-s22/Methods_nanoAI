"""
exp3 3H — restrict the drug signal to the sham->TBI axis.
Build the sham-vs-TBI discriminant (LDA on 66 ratios, standardized) = unit direction w in feature
space. For every drug effect vector v (e.g. resv_TBI - TBI) decompose into ALONG-axis (v.w) and
ORTHOGONAL (|v_perp|); report cos angle = fraction of the drug signal that lies along the disease
axis. Then redo the batch-free within-plate drug test on the 1D axis score only (permutation),
plus DiD and TOST on the axis. This asks: along the sham->TBI direction (ignoring off-axis drug
chemistry), do the drugs move the TBI media?
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
OUT="/mnt/c/Users/riccardo-s/Documents/CNT/TBI/paper/svg/truelog_rerun"; TP=os.environ.get("EXP3_TP","3H")
ORDER=["ch8_3","ch6_5","ch7_5","ch10_2","ch9_4","ch8_4","ch7_6","ch8_6","ch8_7","ch9_5","ch10_3","ch10_5"]
PAIRS=list(combinations(range(12),2)); flo=0.0  # TRUELOG RERUN (was 1e-6)
z=np.load("exp3_%s_peakfeats.npz"%TP)
def ratio(g):
    M=z[f"MX_{g}"]; return np.column_stack([np.log10((M[:,a]+flo)/(M[:,b]+flo)) for a,b in PAIRS])
GR=["sham","tbi","resv_tbi","mg132_tbi","resv_sham","mg132_sham"]
R={g:ratio(g) for g in GR}
colc={"sham":"#2e86c1","resv_sham":"#a9cce3","mg132_sham":"#d2b4de","tbi":"#cb4335","resv_tbi":"#27ae60","mg132_tbi":"#8e44ad"}

# axis from sham vs TBI in standardized ratio space
sc=StandardScaler().fit(np.vstack([R["sham"],R["tbi"]]))
Z={g:sc.transform(R[g]) for g in GR}
lda=LDA().fit(np.vstack([Z["sham"],Z["tbi"]]),np.r_[np.zeros(len(Z["sham"])),np.ones(len(Z["tbi"]))])
w=lda.coef_[0]; flip=1.0 if (Z["tbi"].mean(0)-Z["sham"].mean(0))@w>0 else -1.0
w=flip*w; wn=w/np.linalg.norm(w)                          # unit axis, oriented sham->TBI +
score={g:Z[g]@wn for g in GR}                             # 1D axis coordinate (TBI high)

pr("############ drug signal RESTRICTED to the sham->TBI axis (%s) ############"%TP)
pr("axis 1D scores (mean +/- sd):")
for g in GR: pr("   %-12s % .3f +/- %.3f"%(g,score[g].mean(),score[g].std()))

# ---- decomposition: how much of each drug effect lies ALONG the axis ----
pr("\n=== drug-effect vector decomposition (standardized ratio space) ===")
pr("  %-22s %8s %8s %8s %8s"%("effect vector","|v|","along","|orth|","cos(along)"))
def decomp(base,drug):
    v=Z[drug].mean(0)-Z[base].mean(0); along=v@wn; vn=np.linalg.norm(v); orth=np.sqrt(max(vn**2-along**2,0))
    return vn,along,orth,(along/vn if vn>0 else 0)
for label,base,drug in [("resv on TBI","tbi","resv_tbi"),("MG132 on TBI","tbi","mg132_tbi"),
                        ("resv on sham","sham","resv_sham"),("MG132 on sham","sham","mg132_sham")]:
    vn,al,ot,cs=decomp(base,drug); pr("  %-22s %8.2f %+8.2f %8.2f %+8.2f"%(label,vn,al,ot,cs))
# DiD vector decomposition
def did_vec(dt,ds): return (Z[dt].mean(0)-Z["tbi"].mean(0))-(Z[ds].mean(0)-Z["sham"].mean(0))
for label,dt,ds in [("DiD resveratrol","resv_tbi","resv_sham"),("DiD MG132","mg132_tbi","mg132_sham")]:
    v=did_vec(dt,ds); al=v@wn; vn=np.linalg.norm(v); ot=np.sqrt(max(vn**2-al**2,0))
    pr("  %-22s %8.2f %+8.2f %8.2f %+8.2f"%(label,vn,al,ot,(al/vn if vn>0 else 0)))

# ---- within-plate drug test ON THE AXIS (1D), batch-free, permutation ----
def perm_F_1d(groups,B=10000,seed=0):
    vals=np.concatenate([groups[g] for g in groups]); lab=np.concatenate([[g]*len(groups[g]) for g in groups])
    def F(l):
        gm=vals.mean(); a=len(set(l)); N=len(vals)
        ssb=sum(len(vals[l==g])*(vals[l==g].mean()-gm)**2 for g in set(l))
        ssw=sum(((vals[l==g]-vals[l==g].mean())**2).sum() for g in set(l))
        return (ssb/(a-1))/(ssw/(N-a))
    f0=F(lab); rng=np.random.default_rng(seed); c=sum(F(rng.permutation(lab))>=f0 for _ in range(B))
    return f0,(c+1)/(B+1)
def perm_t_1d(a,b,B=10000,seed=0):
    obs=b.mean()-a.mean(); pooled=np.r_[a,b]; n=len(a); rng=np.random.default_rng(seed)
    c=sum(abs(np.diff([pooled[p][:n].mean(),pooled[p][n:].mean()])[0])>=abs(obs) for p in (rng.permutation(len(pooled)) for _ in range(B)))
    return obs,(c+1)/(B+1)
pr("\n=== within-plate drug test ON AXIS (1D, batch-free, permutation B=10000) ===")
f,p=perm_F_1d({g:score[g] for g in ["tbi","resv_tbi","mg132_tbi"]}); pr("  TBI plate 3-group on axis: F=%.2f p=%.4f"%(f,p))
for a,b in [("tbi","resv_tbi"),("tbi","mg132_tbi")]:
    o,p=perm_t_1d(score[a],score[b]); pr("    %-20s axis-diff=%+.3f p=%.4f"%(a+" vs "+b,o,p))
f,p=perm_F_1d({g:score[g] for g in ["sham","resv_sham","mg132_sham"]}); pr("  SHAM plate 3-group on axis: F=%.2f p=%.4f (ref)"%(f,p))

# ---- DiD on axis + TOST ----
def did_axis(dt,ds):
    a,b,c,d=score["tbi"],score[dt],score["sham"],score[ds]
    did=(b.mean()-a.mean())-(d.mean()-c.mean())
    se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b)+c.var(ddof=1)/len(c)+d.var(ddof=1)/len(d))
    return did,se,2*stats.t.sf(abs(did/se),len(a)+len(b)+len(c)+len(d)-4)
pr("\n=== DiD on axis (neg=rescue) ===")
for label,dt,ds in [("resveratrol","resv_tbi","resv_sham"),("MG132","mg132_tbi","mg132_sham")]:
    did,se,p=did_axis(dt,ds); pr("  %-12s DiD=%+.3f (SE %.3f) p=%.4f"%(label,did,se,p))
gap=score["tbi"].mean()-score["sham"].mean()
def pct_rescue(drug,seed=0):
    rng=np.random.default_rng(seed); b=[]
    for _ in range(5000):
        s=score["sham"][rng.integers(0,len(score["sham"]),len(score["sham"]))].mean()
        t=score["tbi"][rng.integers(0,len(score["tbi"]),len(score["tbi"]))].mean()
        r=score[drug][rng.integers(0,len(score[drug]),len(score[drug]))].mean()
        if abs(t-s)>1e-9: b.append(100*(t-r)/(t-s))
    return np.median(b),np.percentile(b,2.5),np.percentile(b,97.5)
m,lo,hi=pct_rescue("resv_tbi"); pr("  resveratrol %% rescue ALONG AXIS = %.0f%% (95%% CI %.0f..%.0f%%)"%(m,lo,hi))
m,lo,hi=pct_rescue("mg132_tbi"); pr("  MG132        %% move ALONG AXIS = %.0f%% (95%% CI %.0f..%.0f%%)"%(m,lo,hi))

# ---- figure ----
fig,ax=plt.subplots(1,2,figsize=(13,5.5))
order6=["sham","resv_sham","mg132_sham","tbi","resv_tbi","mg132_tbi"]
for xi,g in enumerate(order6):
    v=score[g]; ax[0].scatter(np.full(len(v),xi)+np.random.uniform(-.1,.1,len(v)),v,c=colc[g],s=24,edgecolor="k",lw=.3,zorder=3)
    ax[0].plot([xi-.27,xi+.27],[v.mean()]*2,c="k",lw=2,zorder=4)
ax[0].axhline(score["sham"].mean(),c="#2e86c1",ls="--",lw=1); ax[0].axhline(score["tbi"].mean(),c="#cb4335",ls="--",lw=1)
ax[0].set_xticks(range(6)); ax[0].set_xticklabels(order6,rotation=40,ha="right",fontsize=8)
ax[0].set_ylabel("sham->TBI axis score (1D)"); ax[0].set_title("drug arms on the disease axis only",fontsize=10)
# along vs orth bars
labs=["resv/TBI","MG132/TBI","resv/sham","MG132/sham"]; al=[];ot=[]
for base,drug in [("tbi","resv_tbi"),("tbi","mg132_tbi"),("sham","resv_sham"),("sham","mg132_sham")]:
    vn,a,o,c=decomp(base,drug); al.append(abs(a)); ot.append(o)
xk=np.arange(4)
ax[1].bar(xk,al,0.6,label="|along axis|",color="#2c3e50")
ax[1].bar(xk,ot,0.6,bottom=al,label="|orthogonal|",color="#bdc3c7")
ax[1].set_xticks(xk); ax[1].set_xticklabels(labs,rotation=30,ha="right",fontsize=8)
ax[1].set_ylabel("drug-effect magnitude (std ratio space)"); ax[1].legend(fontsize=8)
ax[1].set_title("most drug signal is OFF the disease axis",fontsize=10)
fig.suptitle("exp3 %s — drug signal restricted to the sham->TBI axis"%TP,fontsize=12)
fig.tight_layout(rect=[0,0,1,0.95]); p=os.path.join(OUT,"exp3_%s_drug_axis.png"%TP); fig.savefig(p,dpi=150); pr("\nsaved",p)
