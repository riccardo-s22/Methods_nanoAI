import os, sys, json, numpy as np, matplotlib
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "shared"))
from repo_paths import results_dir
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# step 2 of 2: reads the JSON written by make_supp_auc_compute.py
OUT=results_dir("01_swcnt_blast_classification","supp_fig2_auc")
d=json.load(open(f"{OUT}/supp_auc_results.json"))
INK,INK2,GRID="#0b0b0b","#52514e","#e5e4e0"; C_SHAM,C_TBI="#2a78d6","#eb6834"
plt.rcParams.update({"font.family":"Arial","font.size":8,"axes.edgecolor":INK2,"axes.labelcolor":INK,
  "xtick.color":INK2,"ytick.color":INK2,"axes.spines.top":False,"axes.spines.right":False,
  "figure.dpi":200,"savefig.bbox":"tight","svg.fonttype":"none","pdf.fonttype":42})
NAMES=["LDA (PCA20)","LogReg L2 (PCA20)","Linear SVM (PCA20)","Elastic-net (PCA20)","Elastic-net (full 2560)"]
SHORT=["LDA","LogReg L2","Linear SVM","Elastic-net\n(PCA20)","Elastic-net\n(2560 feat.)"]

fig=plt.figure(figsize=(7.2,8.3))
gs=fig.add_gridspec(3,2,hspace=0.50,wspace=0.30,height_ratios=[1,1,0.95])

for j,(TAG,lab) in enumerate([("2h","2 h"),("24h","24 h")]):
    ax=fig.add_subplot(gs[0,j]); r=d[TAG]
    ax.plot([0,1],[0,1],lw=0.7,ls=":",color=INK2,zorder=1)
    fpr,tpr=r["roc"][NAMES[0]]
    ax.plot(fpr,tpr,lw=2.0,color="#1b3f8b",solid_capstyle="butt",zorder=4)
    ax.set_xlim(-0.04,1.04); ax.set_ylim(-0.04,1.08)
    ax.set_xlabel("false positive rate"); ax.set_ylabel("true positive rate")
    ax.set_title(f"{'ab'[j]}   leave-one-organoid-out ROC, {lab}",loc="left",fontsize=8.5,fontweight="bold")
    ax.grid(color=GRID,lw=0.5); ax.set_axisbelow(True)
    ax.text(0.97,0.03,"AUC = 1.00 for all five classifiers\n(curves identical, superimposed)\n\n"
            +"\n".join("  · "+n for n in NAMES)
            +f"\n\nn = {r['n_sham']} sham, {r['n_tbi']} blast organoids",
            transform=ax.transAxes,fontsize=6.2,color=INK,ha="right",va="bottom",linespacing=1.35)

ax=fig.add_subplot(gs[1,0])
rng=np.random.default_rng(1); pos=0; xt=[]; xl=[]; gaps={}
for TAG,lab in [("2h","2 h"),("24h","24 h")]:
    s=np.array(d[TAG]["lda_loo_scores"]); y=np.array(d[TAG]["y"])
    gaps[TAG]=float(s[y==1].min()-s[y==0].max())
    for g,cl,nm in [(0,C_SHAM,"sham"),(1,C_TBI,"blast")]:
        v=s[y==g]; ax.scatter(pos+rng.uniform(-.17,.17,len(v)),v,s=11,color=cl,alpha=.8,lw=0,zorder=3)
        ax.plot([pos-.30,pos+.30],[v.mean()]*2,color=INK,lw=1.4,zorder=4)
        xt.append(pos); xl.append(f"{nm}\n{lab}"); pos+=1
    pos+=0.45
ax.axhline(0,color=INK2,lw=.7,ls="--",zorder=2)
ax.set_xticks(xt,xl,fontsize=7); ax.set_ylabel("out-of-fold LDA decision score")
ax.set_title("c   held-out score separation",loc="left",fontsize=8.5,fontweight="bold")
ax.grid(axis="y",color=GRID,lw=.5); ax.set_axisbelow(True)
ax.set_ylim(top=ax.get_ylim()[1]*1.38)
ax.text(.02,.97,"no overlap at either time point\ngap between groups: %.0f (2 h), %.0f (24 h)"%(gaps["2h"],gaps["24h"]),
        transform=ax.transAxes,fontsize=6.4,color=INK2,va="top",linespacing=1.3)

ax=fig.add_subplot(gs[1,1])
bins=np.linspace(0.15,1.02,36)
for TAG,lab,cl in [("2h","2 h","#2a78d6"),("24h","24 h","#7a8699")]:
    ax.hist(d[TAG]["null_well"],bins=bins,histtype="step",lw=1.3,color=cl,label=f"null, {lab} (1,000 shuffles)")
ax.axvline(1.0,color=C_TBI,lw=1.8,zorder=5)
ymax=ax.get_ylim()[1]; ax.set_ylim(0,ymax*1.30)
ax.annotate("observed\nAUC = 1.00\n$P$ = 0.001",xy=(.995,ymax*.52),xytext=(.76,ymax*.80),fontsize=6.8,
            color=INK,ha="center",va="center",linespacing=1.3,
            arrowprops=dict(arrowstyle="->",color=C_TBI,lw=1.0,shrinkA=2,shrinkB=1))
ax.set_xlabel("leave-one-organoid-out AUC"); ax.set_ylabel("permutations")
ax.set_title("d   label-permutation null",loc="left",fontsize=8.5,fontweight="bold")
ax.set_xlim(.15,1.06); ax.grid(axis="y",color=GRID,lw=.5); ax.set_axisbelow(True)
ax.legend(loc="upper left",frameon=False,fontsize=6.3,handlelength=1.6)

# --- e : AUC grid, 5 classifiers x (3 held-out units x 2 timepoints) --------
ax=fig.add_subplot(gs[2,:])
schemes=[("well","organoid"),("row","plate row"),("col","plate column")]
M=np.ones((5,6))
for si,(sk,_) in enumerate(schemes):
    for ti,TAG in enumerate(("2h","24h")):
        for ci,nm in enumerate(NAMES):
            M[ci,si*2+ti]= d[TAG]["auc"][nm] if sk=="well" else d[TAG]["scheme"][sk][nm]["auc"]
ax.imshow(M,cmap="Blues",vmin=0.5,vmax=1.0,aspect="auto")
for i in range(5):
    for j in range(6): ax.text(j,i,"%.2f"%M[i,j],ha="center",va="center",fontsize=7.5,color="white",fontweight="bold")
ax.set_yticks(range(5),SHORT,fontsize=7)
ax.set_xticks(range(6),["2 h","24 h"]*3,fontsize=7.5)
for si,(_,lab) in enumerate(schemes):
    ax.text(si*2+0.5,-0.80,"held out: one %s"%lab,ha="center",va="center",fontsize=7.5,color=INK)
    if si: ax.axvline(si*2-0.5,color="white",lw=3)
ax.set_xlim(-0.5,5.5); ax.set_ylim(4.5,-1.15)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)
ax.set_title("e   AUC is unchanged by what is held out",loc="left",fontsize=8.5,fontweight="bold",pad=16)

for ext in (".svg",".pdf",".png"): plt.savefig(f"{OUT}/supp_auc_blast_vs_sham{ext}")
print("wrote", OUT+"/supp_auc_blast_vs_sham.{svg,pdf,png}")
for t in ("2h","24h"):
    n=np.array(d[t]["null_well"])
    print("  %s: perm p=%.4f | null mean %.3f, 95th %.3f, max %.3f | LDA gap %+.0f"
          % (t,d[t]["perm_p_well"],n.mean(),np.percentile(n,95),n.max(),gaps[t]))
