import glob, numpy as np, pandas as pd, anndata
from scipy.stats import spearmanr, fisher_exact
from statsmodels.stats.multitest import multipletests
# Section path configuration: see ../gf_paths.py. Inputs: the excitatory-neuron
# AnnData below $GF_WORK/h5ad/ (written by 01_dataset_preparation) and the per-arm
# driver lists <arm>_hits_fdr05.tsv written by excit_pathway_heatmap_enrich.py
# (a rerun below <results>/.../pathway_heatmap/enrich, else
# $GF_WORK/results/pathway_heatmap/enrich). Neither is shipped in data/.
# Outputs: <results>/02_geneformer_in_silico_screen/confounders/.
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gf_paths
B=gf_paths.GF_WORK
ENRICH=gf_paths.analysis_input("pathway_heatmap","enrich")
OUTDIR=gf_paths.gf_results_dir("confounders")
YGENES={"UTY","KDM5D","DDX3Y","USP9Y","ZFY","EIF1AY","RPS4Y1","RPS4Y2","NLGN4Y","TXLNGY","PRKY","TTTY14","TTTY15","LINC00278","TMSB4Y","SRY","AMELY","TBL1Y","PCDH11Y","TTTY10","XIST","TSIX","JPX","FTX"}
SEXCHROM_XIST={"XIST","TSIX","JPX","FTX"}
a=anndata.read_h5ad(f"{B}/h5ad/celltype_subsets_paperQC_nGene1000/excitatory_neurons_raw_geneformer_ready.h5ad")
sym=a.var["gene_symbol"].astype(str).values
X=a.X.tocsr()
# detection universe (>=20% in TBI or control)
det=np.maximum(*[np.asarray((X[(a.obs.condition==c).values]>0).mean(0)).ravel() for c in ["TBI","control"]])
univ=set(sym[det>=0.2])
# age correlation across TBI donors (pseudobulk log CPM), donors with >=20 cells
t=a.obs.condition.values=="TBI"
rows,ages=[],[]
for d,g in a.obs[t].groupby("donor",observed=True):
    if len(g)<20: continue
    idx=np.where(a.obs_names.isin(g.index))[0]
    v=np.asarray(X[idx].sum(0)).ravel(); rows.append(np.log1p(v/v.sum()*1e6)); ages.append(float(g.age.iloc[0]))
M=np.vstack(rows); ages=np.array(ages)
keep=det>=0.2
rho=np.full(len(sym),np.nan); p=np.full(len(sym),np.nan)
for j in np.where(keep)[0]:
    r,pp=spearmanr(ages,M[:,j]); rho[j]=r; p[j]=pp
ok=~np.isnan(p); q=np.full(len(sym),np.nan); q[ok]=multipletests(p[ok],method="fdr_bh")[1]
agedf=pd.DataFrame({"gene":sym,"rho":rho,"p":p,"q":q})[keep]
agesig=set(agedf[(agedf.q<0.05)].gene)
agestrong=set(agedf[(agedf.rho.abs()>=0.8)].gene)
print(f"TBI donors used: {len(ages)} ages {sorted(ages)}; universe {len(univ)}; age-assoc FDR<0.05: {len(agesig)}; |rho|>=0.8: {len(agestrong)}")
lists={}
for f in sorted(glob.glob(f"{ENRICH}/*/*_hits_fdr05.tsv")):
    lab=f.split("/")[-2]; lists["pathway:"+lab]=set(pd.read_csv(f,sep="\t").gene_symbol.dropna().astype(str))
out=[]
for k,s in lists.items():
    s&=univ
    y=sorted(s&YGENES); ag=s&agestrong
    # enrichment of strong age genes vs universe
    a_=len(ag); b_=len(s)-a_; c_=len(agestrong&univ)-a_; d_=len(univ)-len(s)-c_
    OR,pf=fisher_exact([[a_,b_],[c_,d_]]) if len(s) else (np.nan,np.nan)
    out.append((k,len(s),",".join(y) or "-",a_,round(100*a_/max(len(s),1),1),round(OR,2),pf))
    print(f"{k:28s} n={len(s):5d} sex-linked={','.join(y) or '-':25s} age|rho|>=0.8: {a_:4d} ({100*a_/max(len(s),1):.1f}%) OR={OR:.2f} p={pf:.2g}")
print("background % strong-age:", round(100*len(agestrong&univ)/len(univ),1))
pd.DataFrame(out,columns=["list","n","sex_linked","n_age_strong","pct","OR","fisher_p"]).to_csv(f"{OUTDIR}/confound_summary.tsv",sep="\t",index=False)
agedf.to_csv(f"{OUTDIR}/age_corr.tsv",sep="\t",index=False)
