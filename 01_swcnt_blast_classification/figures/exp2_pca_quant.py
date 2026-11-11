"""Quantify sham vs blast separation in the shape-spectrum PCA of Fig. 2b / Supp. Fig. 3b (exp2).
PCA exactly as plotted in exp2_analysis2.py: shape vectors (each excitation divided by its sum, 2,560 features),
StandardScaler -> PCA. Reported per read (2 h, 24 h):
  variance explained by PC1/PC2; PC1 sham vs blast: AUC, Mann-Whitney p, Cohen's d
  silhouette with the true labels: in the PC1-PC2 plane and in the full standardized space
  PERMANOVA (Euclidean, full standardized space): pseudo-F, R2, permutation p (B = 9999)
  dispersion check (PERMDISP-like): distance of each well to its group centroid, permutation t-test
  k-means (k = 2, PC1-PC2) without labels: accuracy and ARI vs condition, ARI vs plate row
Seed 0. Run: python exp2_pca_quant.py
(writes results/01_swcnt_blast_classification/pca_quant/exp2_pca_quant_results.txt)
"""
import os, sys
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[v] = "1"
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, ".."))                    # exp2_load.py
sys.path.insert(0, os.path.join(_HERE, "..", "..", "shared"))
from repo_paths import data_path, results_dir
import importlib, numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import roc_auc_score, silhouette_score, adjusted_rand_score as ARI
from scipy.stats import mannwhitneyu
HERE = os.path.dirname(os.path.abspath(__file__)); B = 9999; OUT = []
def pr(s=""): print(s); OUT.append(s)

def permanova(X, y, rng):
    def F(lab):
        m = X.mean(0); sst = ((X - m) ** 2).sum()
        ssw = sum(((X[lab == g] - X[lab == g].mean(0)) ** 2).sum() for g in (0, 1))
        return ((sst - ssw) / 1) / (ssw / (len(X) - 2)), (sst - ssw) / sst
    f0, r2 = F(y); c = sum(F(rng.permutation(y))[0] >= f0 for _ in range(B))
    return f0, r2, (c + 1) / (B + 1)
def disp(X, y, rng):
    d = np.empty(len(y))
    for g in (0, 1): d[y == g] = np.linalg.norm(X[y == g] - X[y == g].mean(0), axis=1)
    obs = d[y == 1].mean() - d[y == 0].mean()
    c = sum(abs((lambda p: d[p == 1].mean() - d[p == 0].mean())(rng.permutation(y))) >= abs(obs) for _ in range(B))
    return d[y == 0].mean(), d[y == 1].mean(), (c + 1) / (B + 1)

for tp in ("2h", "24h"):
    os.environ["EXP2_DIR"] = data_path("01_swcnt_blast_classification", "spectra", tp); os.environ["MIX_INVERTED"] = "1"
    import exp2_load; importlib.reload(exp2_load)
    cube, em = exp2_load.load_all(); rows, cols, grp = exp2_load.labels()
    m = np.isin(grp, ["sham", "TBI"]); y = (grp[m] == "TBI").astype(int); rr = rows[m]; C = cube[m].astype(float)
    S = C.copy()
    for k in range(C.shape[2]): S[:, :, k] /= S[:, :, k].sum(1, keepdims=True)
    Xs = StandardScaler().fit_transform(S.reshape(len(C), -1))
    pca = PCA(2, random_state=0).fit(Xs); Z = pca.transform(Xs); rng = np.random.default_rng(0)
    pr("==== exp2 %s  (sham n=%d, blast n=%d) ====" % (tp, (y == 0).sum(), (y == 1).sum()))
    pr("  variance explained: PC1 %.1f%%, PC2 %.1f%%" % tuple(100 * pca.explained_variance_ratio_))
    for j in (0, 1):
        a, b = Z[y == 0, j], Z[y == 1, j]; auc = roc_auc_score(y, Z[:, j]); auc = max(auc, 1 - auc)
        sd = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
        pr("  PC%d sham vs blast: AUC %.3f, Mann-Whitney p = %.1e, Cohen's d = %+.2f" % (j + 1, auc, mannwhitneyu(a, b).pvalue, (b.mean() - a.mean()) / sd))
    pr("  silhouette (true labels): PC1-PC2 plane %+.3f | full standardized space %+.3f" % (silhouette_score(Z, y), silhouette_score(Xs, y)))
    f, r2, p = permanova(Xs, y, rng)
    pr("  PERMANOVA (full space, Euclidean): pseudo-F = %.2f, R2 = %.3f, p = %.4f (B = %d)" % (f, r2, p, B))
    d0, d1, pd_ = disp(Xs, y, rng)
    pr("  dispersion: mean distance to centroid sham %.1f, blast %.1f, permutation p = %.4f" % (d0, d1, pd_))
    lab = KMeans(2, n_init=25, random_state=0).fit_predict(Z)
    pr("  k-means (k=2, PC1-PC2, no labels): accuracy %.3f, ARI vs condition %+.3f, ARI vs plate row %+.3f"
       % (max((lab == y).mean(), (lab != y).mean()), ARI(y, lab), ARI(rr, lab)))
    pr()
open(os.path.join(results_dir("01_swcnt_blast_classification", "pca_quant"), "exp2_pca_quant_results.txt"), "w").write("\n".join(OUT) + "\n")
