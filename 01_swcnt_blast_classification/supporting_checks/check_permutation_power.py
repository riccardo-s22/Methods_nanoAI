"""Power analysis for the permutation tests behind the blast-vs-sham AUC."""
import os, sys, itertools
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"): os.environ[v]="1"
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, ".."))                    # exp2_load.py
sys.path.insert(0, os.path.join(_HERE, "..", "..", "shared"))
from repo_paths import data_path
import numpy as np, warnings; warnings.filterwarnings("ignore")
from math import comb
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.base import clone
from sklearn.model_selection import LeaveOneGroupOut

print("="*78)
print("A.  MONTE-CARLO RESOLUTION  (well-level test, p = (c+1)/(n+1))")
print("="*78)
print("  %8s %12s %28s" % ("n perms","min p","smallest claim supportable"))
for n in (99,199,999,1999,9999,99999):
    print("  %8d %12.5f %28s" % (n, 1/(n+1), "p < %.3g" % (10**np.ceil(np.log10(1/(n+1))))))
print("\n  MC standard error of p-hat when the true tail prob is p0:  sqrt(p0(1-p0)/n)")
for p0 in (0.05,0.01,0.001):
    print("    p0=%-6g  n=1000 -> SE %.4f   n=10000 -> SE %.5f" % (p0, np.sqrt(p0*(1-p0)/1000), np.sqrt(p0*(1-p0)/10000)))

print("\n"+"="*78)
print("B.  DESIGN CEILING of the EXACT column-level test")
print("="*78)
print("  k columns/group, C(2k,k) assignments, C(2k,k)/2 distinct partitions")
print("  %6s %12s %12s %12s %s" % ("k","assignments","partitions","min p","can reject at..."))
for k in range(2,9):
    N=comb(2*k,k); part=N//2; mp=1/part
    can=[a for a in (0.05,0.01,0.001) if mp<=a]
    print("  %6d %12d %12d %12.4f %s" % (k,N,part,mp, (", ".join("%.3g"%a for a in can)) if can else "NOTHING"))
print("\n  -> your design (k=4): min p = 1/35 = 0.0286.")
print("     Power at alpha=0.05 can reach 1.0; power at alpha=0.01 is EXACTLY 0 for any effect size.")

print("\n"+"="*78)
print("C.  POWER vs EFFECT SIZE, exact 4-vs-4 column test at alpha=0.05")
print("="*78)
print("  At k=4 the only outcome with p<=0.05 is p=1/35, i.e. the 4 blast columns must")
print("  rank strictly above the 4 sham columns. Power = P(perfect separation of 4 vs 4).")
rng=np.random.default_rng(0); NS=200000
print("\n  %8s %10s %10s" % ("d","power","(MC SE)"))
rows=[]
for d in (0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0):
    a=rng.normal(0,1,(NS,4)); b=rng.normal(d,1,(NS,4))
    pw=float((b.min(1)>a.max(1)).mean()); rows.append((d,pw))
    print("  %8.1f %10.3f %10.4f" % (d,pw,np.sqrt(pw*(1-pw)/NS)))
ds=np.array([r[0] for r in rows]); pws=np.array([r[1] for r in rows])
for target in (0.8,0.9):
    lo=np.searchsorted(pws,target)
    if 0<lo<len(ds):
        dd=ds[lo-1]+(target-pws[lo-1])*(ds[lo]-ds[lo-1])/(pws[lo]-pws[lo-1])
        print("  -> d needed for %.0f%% power: %.2f (between-column SD units)" % (target*100, dd))
print("\n  type-I error at d=0 is %.4f, matching the exact 1/35 = %.4f (test is exact, not conservative)"
      % (pws[0], 1/35))

print("\n  Power if the design had k columns per group (d fixed at the value observed below):")
print("\n"+"="*78)
print("D.  OBSERVED column-level effect size in your data")
print("="*78)
for TAG,D in [("2h",data_path("01_swcnt_blast_classification","spectra","2h")),
              ("24h",data_path("01_swcnt_blast_classification","spectra","24h"))]:
    os.environ["EXP2_DIR"]=D
    for m in [k for k in list(sys.modules) if k.startswith("exp2_load")]: del sys.modules[m]
    from exp2_load import load_all, labels
    cube,em=load_all(); rows_,cols,grp=labels()
    mask=np.isin(grp,["sham","TBI"]); y=(grp[mask]=="TBI").astype(int); cc=cols[mask]; C=cube[mask]
    S=C.copy()
    for k in range(C.shape[2]): S[:,:,k]/=S[:,:,k].sum(1,keepdims=True)
    X=S.reshape(len(C),-1)
    mo=make_pipeline(StandardScaler(),PCA(20,random_state=0),LDA())
    s=np.zeros(len(y))
    for tr,te in LeaveOneGroupOut().split(X,y,cc):
        s[te]=clone(mo).fit(X[tr],y[tr]).decision_function(X[te])
    colmean=np.array([s[cc==c].mean() for c in sorted(set(cc))])
    coly=np.array([y[cc==c][0] for c in sorted(set(cc))])
    a,b=colmean[coly==0],colmean[coly==1]
    sp=np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
    d=(b.mean()-a.mean())/sp
    print("  [%s] column means sham %s" % (TAG, np.array2string(a,precision=0)))
    print("       column means blast %s" % np.array2string(b,precision=0))
    print("       between-column pooled SD = %.1f   Cohen's d = %.2f   separation gap = %.1f (%.2f SD)"
          % (sp, d, b.min()-a.max(), (b.min()-a.max())/sp))
    pw=float((rng.normal(0,1,(NS,4)).max(1) < rng.normal(d,1,(NS,4)).min(1)).mean())
    print("       -> power of the exact 4v4 test at this d: %.3f" % pw)
    for k in (3,4,5,6):
        pwk=float((rng.normal(0,1,(NS,k)).max(1) < rng.normal(d,1,(NS,k)).min(1)).mean())
        mp=1/(comb(2*k,k)//2)
        print("          k=%d columns/group: min p=%.4f, power(alpha=.05)=%s"
              % (k, mp, ("%.3f"%pwk) if mp<=0.05 else "0 (cannot reject)"))
