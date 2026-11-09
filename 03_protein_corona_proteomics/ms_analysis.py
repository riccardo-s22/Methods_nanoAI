"""
TBI blast proteomics (Spectronaut PG report) - detected proteome, differential
abundance (blast vs control) and enrichment analysis.

Control: C4R, C5R, C9    Blast: A1, A2, C6R  (C4, C5 excluded: separate batch)
"""
import os, json, time
import numpy as np, pandas as pd
import requests
from scipy import stats, special
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text

BASE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(BASE, "..", "Riccado_Report.csv")
OUT = BASE
FIG = os.path.join(OUT, "figures"); os.makedirs(FIG, exist_ok=True)

CTRL = ["C4R", "C5R", "C9"]
BLAST = ["A1", "A2", "C6R"]
SAMPLES = CTRL + BLAST
# C4/C5 were acquired in a different batch than the "R" re-runs, which cluster with
# all blast samples (Pearson r ~0.976 vs ~0.995). Batch is modelled as a covariate.
BATCH = {"A1": 2, "C4": 1, "C5": 1, "C4R": 2, "C5R": 2, "A2": 2, "C6R": 2, "C9": 2}
FC_CUT, FDR_CUT = 0.0, 0.05            # DAP = BH-adj p < 0.05, no fold-change cutoff

C_CTRL, C_BLAST = "#2a78d6", "#eb6834"
C_UP, C_DOWN, C_NS = "#e34948", "#2a78d6", "#c3c2b7"
INK, INK2 = "#0b0b0b", "#52514e"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight",
                     "svg.fonttype": "none", "pdf.fonttype": 42})   # editable text in SVG/PDF

# ----------------------------------------------------------------- load
raw = pd.read_csv(RAW)
qcols = {c: c.split("] ")[1].split(".raw")[0] for c in raw.columns if "Quantity" in c}
raw = raw.rename(columns=qcols)
ALL_RUNS = list(qcols.values())

def pick_canonical(row):
    """Representative accession: first Swiss-Prot entry in the group (mnemonic name)."""
    accs = row["PG.ProteinGroups"].split(";"); names = row["PG.ProteinNames"].split(";")
    descs = row["PG.ProteinDescriptions"].split(";")
    for i, (a, n) in enumerate(zip(accs, names)):
        if n.split("_")[0] != a.split("-")[0]:
            return a, n, descs[min(i, len(descs) - 1)]
    return accs[0], names[0], descs[0]

rep = raw.apply(pick_canonical, axis=1, result_type="expand")
raw["Accession"], raw["EntryName"], raw["Description"] = rep[0], rep[1], rep[2]

# UniProt accession -> gene symbol via g:Profiler g:Convert
def gconvert(ids):
    r = requests.post("https://biit.cs.ut.ee/gprofiler/api/convert/convert/",
                      json={"organism": "hsapiens", "query": list(ids), "target": "ENSG"}, timeout=300)
    r.raise_for_status()
    hits = {}
    for e in r.json()["result"]:
        if e["name"] not in (None, "None", "N/A"):
            hits.setdefault(e["incoming"], []).append((e["name"], e["converted"]))
    return hits

def uniprot_genes(accs):
    """Primary gene symbol per accession from UniProt (authoritative for the protein)."""
    out = {}; accs = [a.split("-")[0] for a in accs]
    for i in range(0, len(accs), 80):
        chunk = accs[i:i + 80]
        q = " OR ".join(f"accession:{a}" for a in chunk)
        for attempt in range(3):
            try:
                r = requests.get("https://rest.uniprot.org/uniprotkb/search", timeout=120,
                                 params={"query": q, "fields": "accession,gene_primary", "format": "tsv", "size": 500})
                r.raise_for_status(); break
            except Exception as ex:
                print("uniprot retry", str(ex)[:80]); time.sleep(5)
        else:
            raise RuntimeError("UniProt lookup failed")
        for line in r.text.strip().split("\n")[1:]:
            p = line.split("\t")
            if len(p) > 1 and p[1]: out[p[0]] = p[1].split(";")[0].strip()
    return out

accs = list(raw["Accession"].unique())
up_sym = uniprot_genes(accs)
hits = gconvert(accs)
sym, ensg = {}, {}
for a in accs:
    g = up_sym.get(a.split("-")[0]); h = hits.get(a, [])
    match = [x for x in h if x[0] == g]          # Ensembl gene whose symbol matches UniProt's
    best = match[0] if match else (h[0] if h else (None, None))
    sym[a] = g or best[0]; ensg[a] = best[1]
raw["Gene"] = raw["Accession"].map(sym); raw["ENSG"] = raw["Accession"].map(ensg)
raw["Gene"] = raw["Gene"].fillna(raw["EntryName"].str.replace("_HUMAN", ""))

# ----------------------------------------------------------------- 1. detected proteome
det = raw[["Accession", "Gene", "ENSG", "EntryName", "Description", "PG.ProteinGroups", "PG.MolecularWeight"]].copy()
det["n_runs_detected_(of_27)"] = raw[ALL_RUNS].notna().sum(axis=1)
det["n_detected_selected_(of_6)"] = raw[SAMPLES].notna().sum(axis=1)
det["n_ctrl_(of_3)"] = raw[CTRL].notna().sum(axis=1)
det["n_blast_(of_3)"] = raw[BLAST].notna().sum(axis=1)
det["mean_log2_intensity_selected"] = np.log2(raw[SAMPLES]).mean(axis=1)
for s in ALL_RUNS:
    det[s] = raw[s]
det = det.sort_values("mean_log2_intensity_selected", ascending=False)

qc = pd.DataFrame({"run": ALL_RUNS, "proteins_quantified": [raw[s].notna().sum() for s in ALL_RUNS],
                   "median_intensity": [raw[s].median() for s in ALL_RUNS]})
qc["used_in_DA"] = qc["run"].isin(SAMPLES)
qc["group"] = qc["run"].map({**{c: "Control" for c in CTRL}, **{b: "Blast" for b in BLAST}}).fillna("")

# ----------------------------------------------------------------- 2. differential abundance
lg = np.log2(raw[SAMPLES])
lg = lg - lg.median() + lg.median().median()               # median normalisation
keep = (raw[CTRL].notna().sum(axis=1) >= 2) & (raw[BLAST].notna().sum(axis=1) >= 2)
print(f"proteins passing filter (>=2/3 ctrl & >=2/3 blast): {keep.sum()}")

def limma_fit(Y, X, coef):
    """Per-protein OLS with missing values (as limma lmFit) + empirical Bayes moderation (eBayes)."""
    n = Y.shape[0]
    beta = np.full(n, np.nan); s2 = np.full(n, np.nan); df = np.zeros(n); vu = np.full(n, np.nan)
    for i in range(n):
        y = Y[i]; m = ~np.isnan(y); Xi = X[m]
        Xi = Xi[:, np.abs(Xi).sum(0) > 0]
        rk = np.linalg.matrix_rank(Xi)
        if rk < X.shape[1] or m.sum() - rk < 1:            # coefficient not estimable
            continue
        XtXi = np.linalg.inv(Xi.T @ Xi); b = XtXi @ Xi.T @ y[m]
        r = y[m] - Xi @ b
        beta[i] = b[coef]; df[i] = m.sum() - rk; s2[i] = r @ r / df[i]; vu[i] = XtXi[coef, coef]
    ok = df > 0
    # fitFDist (Smyth 2004): moment estimates of d0, s0^2 on log scale
    z = np.log(s2[ok]); d = df[ok]
    e = z - special.digamma(d / 2) + np.log(d / 2)
    emean = e.mean(); evar = ((e - emean) ** 2).sum() / (len(e) - 1) - special.polygamma(1, d / 2).mean()
    if evar > 0:
        # invert trigamma by Newton
        yv = evar; x = 0.5 + 1 / yv
        for _ in range(50):
            tri = special.polygamma(1, x); dif = tri * (1 - tri / yv) / special.polygamma(2, x)
            x += dif
            if -dif / x < 1e-8: break
        d0 = 2 * x; s02 = np.exp(emean + special.digamma(d0 / 2) - np.log(d0 / 2))
    else:
        d0 = np.inf; s02 = np.exp(emean)
    s2p = (d0 * s02 + df * s2) / (d0 + df) if np.isfinite(d0) else np.full(n, s02)
    t = beta / np.sqrt(s2p * vu)
    dft = df + d0 if np.isfinite(d0) else np.full(n, 1e6)
    p = 2 * stats.t.sf(np.abs(t), dft)
    return beta, t, p, df, d0, s02

Y = lg[keep].to_numpy()
grp = np.array([1.0 if s in BLAST else 0.0 for s in SAMPLES])
bat = np.array([1.0 if BATCH[s] == 1 else 0.0 for s in SAMPLES])

# Primary: ~ group (+ batch when C4/C5 are included)
X1 = np.column_stack([np.ones(len(SAMPLES)), grp, bat]) if bat.any() else np.column_stack([np.ones(len(SAMPLES)), grp])
b1, t1, p1, df1, d01, s01 = limma_fit(Y, X1, 1)
# Unadjusted (as literally requested, no batch term) - for transparency
X0 = np.column_stack([np.ones(len(SAMPLES)), grp])
b0, t0, p0, *_ = limma_fit(Y, X0, 1)
# Sensitivity: batch-matched only (C4R, C5R, C9 vs A1, A2, C6R)
bm = np.array([BATCH[s] == 2 for s in SAMPLES])
Ybm = Y[:, bm]; Xbm = np.column_stack([np.ones(bm.sum()), grp[bm]])
b2, t2, p2, *_ = limma_fit(Ybm, Xbm, 1)
print(f"eBayes prior: d0={d01:.2f}, s0^2={s01:.4f}")

def bh(p):
    q = np.full_like(p, np.nan); ok = ~np.isnan(p)
    q[ok] = multipletests(p[ok], method="fdr_bh")[1]; return q

da = det.loc[keep[keep].index, ["Accession", "Gene", "ENSG", "EntryName", "Description", "PG.ProteinGroups"]].copy()
da = da.loc[lg[keep].index]
da["log2FC_blast_vs_ctrl"] = b1; da["t_moderated"] = t1; da["P.Value"] = p1; da["adj.P.Val_BH"] = bh(p1)
da["mean_log2_ctrl"] = lg.loc[keep, CTRL].mean(axis=1); da["mean_log2_blast"] = lg.loc[keep, BLAST].mean(axis=1)
da["n_ctrl"] = raw.loc[keep, CTRL].notna().sum(axis=1); da["n_blast"] = raw.loc[keep, BLAST].notna().sum(axis=1)
da["log2FC_noBatchTerm"] = b0; da["adj.P_noBatchTerm"] = bh(p0)
da["log2FC_batchMatched"] = b2; da["adj.P_batchMatched"] = bh(p2)
da["DAP"] = np.where((da["adj.P.Val_BH"] < FDR_CUT) & (da["log2FC_blast_vs_ctrl"] >= FC_CUT), "Up",
              np.where((da["adj.P.Val_BH"] < FDR_CUT) & (da["log2FC_blast_vs_ctrl"] <= -FC_CUT), "Down", "ns"))
da["confirmed_batchMatched"] = (da["DAP"] != "ns") & (da["adj.P_batchMatched"] < FDR_CUT) & \
                               (np.sign(da["log2FC_batchMatched"]) == np.sign(da["log2FC_blast_vs_ctrl"]))
for s in SAMPLES:
    da[f"log2norm_{s}"] = lg.loc[keep, s]
da = da.sort_values("P.Value")

# presence/absence (not testable): all ctrl present & no blast, or vice versa
pa = det[["Accession", "Gene", "Description"]].copy()
pa["n_ctrl"] = raw[CTRL].notna().sum(axis=1); pa["n_blast"] = raw[BLAST].notna().sum(axis=1)
pa["mean_log2_ctrl"] = lg[CTRL].mean(axis=1); pa["mean_log2_blast"] = lg[BLAST].mean(axis=1)
onoff = pa[((pa.n_ctrl == len(CTRL)) & (pa.n_blast == 0)) | ((pa.n_ctrl == 0) & (pa.n_blast == len(BLAST)))].copy()
onoff["pattern"] = np.where(onoff.n_blast == 0, "ctrl-only", "blast-only")

summ = {"proteins_in_report": len(raw), "detected_in_any_selected_sample": int((raw[SAMPLES].notna().any(axis=1)).sum()),
        "tested": int(keep.sum()), "Up": int((da.DAP == "Up").sum()), "Down": int((da.DAP == "Down").sum()),
        "Up_conf_batchMatched": int(((da.DAP == "Up") & da.confirmed_batchMatched).sum()),
        "Down_conf_batchMatched": int(((da.DAP == "Down") & da.confirmed_batchMatched).sum()),
        "DAP_noBatchTerm_(FDR<0.05)": int(((da["adj.P_noBatchTerm"] < FDR_CUT) & (da.log2FC_noBatchTerm.abs() >= FC_CUT)).sum()),
        "DAP_batchMatched_(FDR<0.05)": int(((da["adj.P_batchMatched"] < FDR_CUT) & (da["log2FC_batchMatched"].abs() >= FC_CUT)).sum()),
        "nominal_p<0.05": int((da["P.Value"] < 0.05).sum()),
        "ctrl_only": int((onoff.pattern == "ctrl-only").sum()), "blast_only": int((onoff.pattern == "blast-only").sum()),
        "eBayes_d0": round(float(d01), 2)}
print(json.dumps(summ, indent=1))

# ----------------------------------------------------------------- enrichment (g:Profiler)
SOURCES = ["GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC", "WP", "CORUM", "HPA"]
def gost(query, background=None, name=""):
    query = pd.Series(list(query)).dropna().unique()
    if background is not None:
        background = pd.Series(list(background)).dropna().unique()
    if len(query) == 0:
        return pd.DataFrame()
    body = {"organism": "hsapiens", "query": list(query), "sources": SOURCES, "user_threshold": 0.05,
            "significance_threshold_method": "g_SCS", "no_evidences": False, "all_results": False}
    if background is not None:
        body["background"] = list(background); body["domain_scope"] = "custom"
    for attempt in range(3):
        try:
            r = requests.post("https://biit.cs.ut.ee/gprofiler/api/gost/profile/", json=body, timeout=600)
            r.raise_for_status(); break
        except Exception as ex:
            print("retry", name, ex); time.sleep(10)
    res = r.json()
    df = pd.DataFrame(res["result"])
    if df.empty:
        return df
    inv = dict(zip(raw["ENSG"], raw["Gene"]))
    qens = res["meta"]["genes_metadata"]["query"]["query_1"]["ensgs"]
    df["genes"] = df["intersections"].apply(lambda ix: ";".join(sorted({inv.get(qens[k], qens[k]) for k, ev in enumerate(ix) if ev})))
    df["gene_ratio"] = df["intersection_size"] / df["query_size"]
    df["fold_enrichment"] = df["gene_ratio"] / (df["term_size"] / df["effective_domain_size"])
    df = df[["source", "native", "name", "p_value", "term_size", "query_size", "intersection_size",
             "effective_domain_size", "fold_enrichment", "genes"]].rename(columns={"native": "term_id", "p_value": "adj_p_gSCS"})
    return df.sort_values(["source", "adj_p_gSCS"])

# 1b. whole detected proteome vs whole genome; also high-abundance (top 10%) vs detected
enr_all = gost(raw["ENSG"], name="all")
top = det.head(int(len(det) * 0.10))["ENSG"]
enr_top = gost(top, background=raw["ENSG"], name="top10")
# 3. DAP enrichment, background = tested proteins
bg = da["ENSG"]
up = da.loc[da.DAP == "Up", "ENSG"]; dn = da.loc[da.DAP == "Down", "ENSG"]
enr_up = gost(up, bg, "up"); enr_dn = gost(dn, bg, "down"); enr_dap = gost(pd.concat([up, dn]), bg, "all DAP")
# (relaxed list in case the strict lists are small)
rel = da[(da["P.Value"] < 0.05) & (da.log2FC_blast_vs_ctrl.abs() >= 0.58)]
enr_rel_up = gost(rel.loc[rel.log2FC_blast_vs_ctrl > 0, "ENSG"], bg, "rel up")
enr_rel_dn = gost(rel.loc[rel.log2FC_blast_vs_ctrl < 0, "ENSG"], bg, "rel dn")

# 3b. pre-ranked GSEA on moderated t (uses every tested protein, threshold-free)
import gseapy as gp
rnk = da.groupby("Gene")["t_moderated"].mean().sort_values(ascending=False)
gsea_tabs = {}
for lib in ["GO_Biological_Process_2025", "Reactome_Pathways_2024", "KEGG_2021_Human"]:
    try:
        pre = gp.prerank(rnk=rnk, gene_sets=lib, min_size=10, max_size=500, permutation_num=1000,
                         seed=1, threads=4, outdir=None, verbose=False)
        t = pre.res2d.copy(); t["library"] = lib
        gsea_tabs[lib] = t
    except Exception as ex:
        print("GSEA failed", lib, ex)
gsea = pd.concat(gsea_tabs.values()) if gsea_tabs else pd.DataFrame()
if not gsea.empty:
    for c in ["NES", "NOM p-val", "FDR q-val", "FWER p-val"]:
        gsea[c] = pd.to_numeric(gsea[c], errors="coerce")
    gsea = gsea.sort_values("FDR q-val")

# ----------------------------------------------------------------- figures
def savefig(name):
    for ext in (".png", ".pdf", ".svg"):           # retry: Windows-mounted drive can transiently lock files
        for attempt in range(5):
            try:
                plt.savefig(os.path.join(FIG, name + ext)); break
            except OSError:
                time.sleep(2)
    plt.close()

# QC: proteins per run
fig, ax = plt.subplots(figsize=(9, 3))
cols = [C_CTRL if r in CTRL else C_BLAST if r in BLAST else C_NS for r in qc.run]
ax.bar(qc.run, qc.proteins_quantified, color=cols, width=0.7)
ax.set_ylabel("Protein groups quantified"); ax.set_title("Proteins quantified per run (colored = used in DA: blue control, orange blast)", loc="left", fontsize=9)
plt.xticks(rotation=90); ax.grid(axis="y", color="#e5e4e0", lw=0.6); ax.set_axisbelow(True)
savefig("01_QC_proteins_per_run")

# correlation heatmap
corr = lg.dropna().corr()
fig, ax = plt.subplots(figsize=(4.6, 4))
im = ax.imshow(corr, cmap="Blues", vmin=corr.values.min(), vmax=1)
ax.set_xticks(range(len(SAMPLES)), SAMPLES, rotation=90); ax.set_yticks(range(len(SAMPLES)), SAMPLES)
for i in range(len(SAMPLES)):
    for j in range(len(SAMPLES)):
        ax.text(j, i, f"{corr.iat[i, j]:.3f}", ha="center", va="center", fontsize=6,
                color="white" if corr.iat[i, j] > (corr.values.min() + 1) / 2 else INK)
plt.colorbar(im, fraction=0.046, label="Pearson r (log2)"); ax.set_title("Sample correlation", loc="left")
savefig("02_sample_correlation")

# PCA (complete cases)
M = lg.dropna(); Mc = (M.T - M.T.mean()).to_numpy()
U, S, Vt = np.linalg.svd(Mc, full_matrices=False); pcs = U * S; ve = S ** 2 / (S ** 2).sum()
fig, ax = plt.subplots(figsize=(4.5, 4))
for i, s in enumerate(SAMPLES):
    c = C_BLAST if s in BLAST else C_CTRL
    ax.scatter(pcs[i, 0], pcs[i, 1], s=70, color=c, edgecolor="white", lw=2, marker="o" if BATCH[s] == 2 else "s", zorder=3)
    ax.annotate(s, (pcs[i, 0], pcs[i, 1]), xytext=(6, 4), textcoords="offset points", fontsize=8, color=INK)
ax.set_xlabel(f"PC1 ({ve[0]*100:.1f}%)"); ax.set_ylabel(f"PC2 ({ve[1]*100:.1f}%)")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=C_CTRL, label="Control"),
                   Line2D([], [], marker="o", ls="", color=C_BLAST, label="Blast"),
                   ] + ([Line2D([], [], marker="s", ls="", color=INK2, label="Batch 1 (C4, C5)")] if "C4" in SAMPLES else []), frameon=False, fontsize=7)
ax.set_title(f"PCA, {len(M)} complete proteins", loc="left"); ax.grid(color="#e5e4e0", lw=0.6); ax.set_axisbelow(True)
savefig("03_PCA")

# volcano
fig, ax = plt.subplots(figsize=(6, 5))
y = -np.log10(da["P.Value"])
# dot colour = nominal p < 0.05 & |log2FC| > 1 (red up / blue down); dotted line = BH FDR 0.05 threshold
nom = np.where((da["P.Value"] < 0.05) & (da.log2FC_blast_vs_ctrl.abs() > 1), np.where(da.log2FC_blast_vs_ctrl > 0, "Up", "Down"), "ns")
for lab, c in [("ns", C_NS), ("Down", C_DOWN), ("Up", C_UP)]:
    m = nom == lab
    ax.scatter(da.loc[m, "log2FC_blast_vs_ctrl"], y[m], s=8 if lab == "ns" else 12, color=c, alpha=0.6 if lab == "ns" else 0.85,
               lw=0, label=f"{lab}: p<0.05, |log2FC|>1 ({m.sum()})" if lab != "ns" else f"ns ({m.sum()})")
ax.axvline(1, color=INK2, lw=0.6, ls="--"); ax.axvline(-1, color=INK2, lw=0.6, ls="--")
if False:
    ax.axvline(FC_CUT, color=INK2, lw=0.6, ls="--"); ax.axvline(-FC_CUT, color=INK2, lw=0.6, ls="--")
sig = da[da["adj.P.Val_BH"] < FDR_CUT]
if len(sig):
    thr = sig["P.Value"].max(); ax.axhline(-np.log10(thr), color=INK2, lw=0.6, ls=":")
lab = da[da.DAP != "ns"].nsmallest(25, "P.Value")
if len(lab) < 10:   # nothing passes FDR: label the strongest nominal hits instead
    lab = pd.concat([lab, da.nsmallest(12, "P.Value"), da[(da["P.Value"] < 0.05) & (da.log2FC_blast_vs_ctrl.abs() >= 0.8)]]).drop_duplicates("Gene")
texts = [ax.text(r.log2FC_blast_vs_ctrl, -np.log10(r["P.Value"]), r.Gene, fontsize=6.5, color=INK) for _, r in lab.iterrows()]
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.4))
ax.set_xlabel("log2 fold change (Blast / Control)"); ax.set_ylabel("-log10 p (moderated t)")
ax.set_title("Blast (A1, A2, C6R) vs Control (C4R, C5R, C9)\ncolored: p < 0.05 & |log2FC| > 1; dotted line: BH FDR = 0.05", loc="left", fontsize=9)
ax.legend(frameon=False, fontsize=7); savefig("04_volcano")

# heatmap of DAPs (z-score)
dap = da[da.DAP != "ns"]
hm_title = f"Top {min(60, len(dap))} DAPs"
if len(dap) < 10:   # too few DAPs: show the top 50 proteins by nominal p instead
    dap = da.nsmallest(50, "P.Value"); hm_title = "Top 50 proteins by p-value (not all FDR-significant)"
if len(dap) >= 2:
    show = dap.nsmallest(60, "P.Value")
    Z = show[[f"log2norm_{s}" for s in SAMPLES]].copy(); Z.columns = SAMPLES
    # remove batch shift for display: centre C4/C5 onto C4R/C5R level per protein
    if {"C4", "C5"} <= set(SAMPLES):
        sh = Z[["C4", "C5"]].mean(axis=1) - Z[["C4R", "C5R"]].mean(axis=1)
        Z[["C4", "C5"]] = Z[["C4", "C5"]].sub(sh, axis=0)
    Z = Z.sub(Z.mean(axis=1), axis=0).div(Z.std(axis=1), axis=0)
    Z = Z.loc[show.sort_values("log2FC_blast_vs_ctrl").index]
    fig, ax = plt.subplots(figsize=(3.8, 0.13 * len(Z) + 1.2))
    im = ax.imshow(Z, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(len(SAMPLES)), SAMPLES, rotation=90); ax.set_yticks(range(len(Z)), show.loc[Z.index, "Gene"], fontsize=6)
    ax.axvline(len(CTRL) - 0.5, color="white", lw=2)
    plt.colorbar(im, fraction=0.05, label="row z-score")
    ax.set_title(hm_title, loc="left", fontsize=8); savefig("05_DAP_heatmap")

def dotplot(df, title, fname, n=20):
    if df is None or df.empty: return
    d = df.copy(); d = d[d.source != "HPA"] if (d.source != "HPA").any() else d
    d = d.nsmallest(n, "adj_p_gSCS").iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, 0.28 * len(d) + 1))
    sc = ax.scatter(-np.log10(d.adj_p_gSCS), range(len(d)), s=d.intersection_size.clip(upper=400) / d.intersection_size.max() * 180 + 15,
                    c=d.fold_enrichment, cmap="Blues", vmin=0, edgecolor=INK2, lw=0.4)
    lbl = [f"{s} · {nm[:60]}" for s, nm in zip(d.source, d.name)]
    ax.set_yticks(range(len(d)), lbl, fontsize=7); ax.set_xlabel("-log10 adjusted p (g:SCS)")
    plt.colorbar(sc, fraction=0.04, label="fold enrichment"); ax.set_title(title, loc="left", fontsize=9)
    ax.grid(axis="x", color="#e5e4e0", lw=0.6); ax.set_axisbelow(True); savefig(fname)

dotplot(enr_all, "Detected proteome vs human genome (top 20 terms)", "06_enrich_detected_proteome")
dotplot(enr_top, "Top-10% most abundant vs detected proteome", "07_enrich_top10pct_abundant")
dotplot(enr_up, "Up in blast (FDR<0.05), background = tested", "08_enrich_DAP_up")
dotplot(enr_dn, "Down in blast (FDR<0.05), background = tested", "09_enrich_DAP_down")
dotplot(enr_dap, "All DAPs, background = tested", "10_enrich_DAP_all")
dotplot(enr_rel_up, "Exploratory: nominal p<0.05 & log2FC≥0.58 (higher in blast)", "10b_enrich_relaxed_up")
dotplot(enr_rel_dn, "Exploratory: nominal p<0.05 & log2FC≤-0.58 (lower in blast)", "10c_enrich_relaxed_down")

if not gsea.empty:
    g = gsea[gsea["FDR q-val"] < 0.25].copy()
    g = pd.concat([g[g.NES > 0].nsmallest(12, "FDR q-val"), g[g.NES < 0].nsmallest(12, "FDR q-val")]).sort_values("NES")
    if len(g):
        fig, ax = plt.subplots(figsize=(7, 0.28 * len(g) + 1))
        ax.barh(range(len(g)), g.NES, color=[C_UP if v > 0 else C_DOWN for v in g.NES], height=0.7)
        ax.set_yticks(range(len(g)), [t[:70] for t in g.Term], fontsize=7); ax.axvline(0, color=INK2, lw=0.6)
        ax.set_xlabel("NES (positive = higher in blast)")
        ax.set_title("Pre-ranked GSEA (moderated t), FDR q < 0.25", loc="left", fontsize=9)
        ax.grid(axis="x", color="#e5e4e0", lw=0.6); ax.set_axisbelow(True); savefig("11_GSEA")

# ----------------------------------------------------------------- write workbook
xl = os.path.join(OUT, "TBI_blast_proteomics_results.xlsx")
with pd.ExcelWriter(xl, engine="openpyxl") as w:
    pd.DataFrame([summ]).T.rename(columns={0: "value"}).to_excel(w, sheet_name="Summary")
    qc.to_excel(w, sheet_name="QC_runs", index=False)
    det.to_excel(w, sheet_name="1_Detected_proteins", index=False)
    enr_all.to_excel(w, sheet_name="1_Enrich_detected_vs_genome", index=False)
    enr_top.to_excel(w, sheet_name="1_Enrich_top10pct_abundant", index=False)
    da.to_excel(w, sheet_name="2_DA_blast_vs_ctrl", index=False)
    da[da.DAP != "ns"].sort_values("log2FC_blast_vs_ctrl", ascending=False).to_excel(w, sheet_name="2_DAPs", index=False)
    onoff.to_excel(w, sheet_name="2_On-Off_proteins", index=False)
    enr_up.to_excel(w, sheet_name="3_Enrich_DAP_up", index=False)
    enr_dn.to_excel(w, sheet_name="3_Enrich_DAP_down", index=False)
    enr_dap.to_excel(w, sheet_name="3_Enrich_DAP_all", index=False)
    enr_rel_up.to_excel(w, sheet_name="3_Enrich_relaxed_up", index=False)
    enr_rel_dn.to_excel(w, sheet_name="3_Enrich_relaxed_down", index=False)
    if not gsea.empty:
        gsea.to_excel(w, sheet_name="3_GSEA_preranked", index=False)
json.dump(summ, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
for nm, df in [("up", enr_up), ("down", enr_dn), ("all_dap", enr_dap), ("rel_up", enr_rel_up), ("rel_dn", enr_rel_dn),
               ("detected", enr_all), ("top10", enr_top)]:
    print(f"\n=== {nm}: {len(df)} terms");
    if len(df): print(df.nsmallest(12, "adj_p_gSCS")[["source", "name", "adj_p_gSCS", "intersection_size", "fold_enrichment"]].to_string(index=False))
if not gsea.empty:
    print("\n=== GSEA FDR<0.05:", (gsea["FDR q-val"] < 0.05).sum())
    print(gsea[gsea["FDR q-val"] < 0.05].head(30)[["library", "Term", "NES", "FDR q-val"]].to_string(index=False))
print("\nTop DAPs:"); print(da[da.DAP != "ns"].head(40)[["Gene", "log2FC_blast_vs_ctrl", "adj.P.Val_BH", "n_ctrl", "n_blast", "confirmed_batchMatched"]].to_string(index=False))
