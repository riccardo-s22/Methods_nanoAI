"""Can simple descriptive statistics replace the classifier?
Test unsupervised / label-free summaries: PC1-PC3 (no labels used), emission centroid,
total intensity. Then the dilution trend on the same axes."""
import os, sys
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"): os.environ[v]="1"
os.environ.setdefault("MIX_INVERTED","1")
sys.path.insert(0,"/mnt/c/Users/riccardo-s/Documents/CNT/TBI/exp2")
import numpy as np, warnings; warnings.filterwarnings("ignore")
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
from scipy import stats

FRAC={"sham":0.0,"mix25":0.25,"mix50":0.5,"mix75":0.75,"TBI":1.0}
for TAG,D in [("2h","/mnt/c/Users/riccardo-s/Documents/CNT/TBI/exp2/2h"),
              ("24h","/mnt/c/Users/riccardo-s/Documents/CNT/TBI/exp2/24h")]:
    os.environ["EXP2_DIR"]=D
    for m in [k for k in list(sys.modules) if k.startswith("exp2_load")]: del sys.modules[m]
    from exp2_load import load_all, labels, EXC
    cube,em=load_all(); rows,cols,grp=labels()
    keep=np.isin(grp,["sham","TBI","mix25","mix50","mix75"]) & (np.arange(96)!=0)
    g=grp[keep]; cc=cols[keep]; C=cube[keep]
    S=C.copy()
    for k in range(C.shape[2]): S[:,:,k]/=S[:,:,k].sum(1,keepdims=True)
    X=S.reshape(len(C),-1)
    st=np.isin(g,["sham","TBI"]); y=(g[st]=="TBI").astype(int)
    print("\n"+"="*76); print(f"{TAG}   n={len(g)} wells ({st.sum()} sham/TBI + {len(g)-st.sum()} mixes)"); print("="*76)

    # ---- unsupervised PCA (labels never used) ----
    Z=PCA(5,random_state=0).fit_transform(StandardScaler().fit_transform(X))
    print("  UNSUPERVISED PCA (no labels used in fitting)")
    print("   %-6s %9s %9s %10s %10s" % ("axis","AUC","Cohen d","MWU p","t-test p"))
    for k in range(3):
        v=Z[st,k]; a,b=v[y==0],v[y==1]
        sp=np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
        auc=roc_auc_score(y,v); auc=max(auc,1-auc)
        print("   PC%-4d %9.4f %9.2f %10.2g %10.2g" % (k+1,auc,abs(b.mean()-a.mean())/sp,
              stats.mannwhitneyu(a,b).pvalue, stats.ttest_ind(a,b).pvalue))
    # dilution trend on PC1
    for k in range(2):
        fr=np.array([FRAC[x] for x in g]); v=Z[:,k]
        if np.corrcoef(v[st],y)[0,1]<0: v=-v
        rho,pr=stats.spearmanr(fr,v)
        gr=[v[g==n] for n in ["sham","mix25","mix50","mix75","TBI"]]
        F,pa=stats.f_oneway(*gr)
        print("   PC%d dilution series: Spearman rho=%+.3f (p=%.2g) | one-way ANOVA F=%.1f p=%.2g | means %s"
              % (k+1,rho,pr,F,pa,np.array2string(np.array([x.mean() for x in gr]),precision=1)))

    # ---- simple physical descriptors ----
    print("\n  SIMPLE DESCRIPTORS (sham vs TBI)")
    print("   %-26s %9s %9s %10s" % ("descriptor","AUC","Cohen d","MWU p"))
    desc={}
    tot=C.sum((1,2)); desc["total intensity"]=tot
    for k,ex in enumerate(EXC):
        sp_=C[:,:,k]; desc["centroid ex%d"%ex]=(sp_*em).sum(1)/sp_.sum(1)
    allc=np.mean([desc["centroid ex%d"%e] for e in EXC],axis=0); desc["mean centroid (5 exc)"]=allc
    for nm,v in desc.items():
        a,b=v[st][y==0],v[st][y==1]
        sp=np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
        auc=roc_auc_score(y,v[st]); auc=max(auc,1-auc)
        print("   %-26s %9.4f %9.2f %10.2g" % (nm,auc,abs(b.mean()-a.mean())/sp,stats.mannwhitneyu(a,b).pvalue))
    fr=np.array([FRAC[x] for x in g])
    rho,pr=stats.spearmanr(fr,allc)
    gr=[allc[g==n] for n in ["sham","mix25","mix50","mix75","TBI"]]
    F,pa=stats.f_oneway(*gr)
    print("   mean-centroid dilution series: Spearman rho=%+.3f (p=%.2g) | ANOVA F=%.1f p=%.2g | means %s"
          % (rho,pr,F,pa,np.array2string(np.array([x.mean() for x in gr]),precision=2)))
