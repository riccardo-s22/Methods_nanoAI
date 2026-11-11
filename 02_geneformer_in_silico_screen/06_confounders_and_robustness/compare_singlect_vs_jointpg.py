#!/usr/bin/env python
"""Compare the single-cell-type ZERO-SHOT ISP (pretrained Geneformer-V2, per-CT
control/TBI centroids) against the JointPG MULTI-TASK fine-tuned ISP, on the
SAME Garza dataset and the SAME TBI->control delete arm, for oligo / astro /
microglia.

Two questions:
 (1) Concordance: do the two models agree on the driver ranking per cell type?
     hit-set overlap (paper-threshold), top-N overlap, Spearman on shared genes.
 (2) Does the JointPG mechano/ECM/Polycomb signal also appear in the zero-shot
     single-CT hits? mechano-panel + ECM-Organization Fisher; consensus-PcG
     Fisher (same sets used in the three-layer analysis).

Locations (see ../gf_paths.py): zero-shot tables below $GF_WORK/results/, MT
tables below $GF_JOINTPG/results/isp/, the mechano gene list below $GF_MECHANO;
none of them is shipped in data/. Output:
<results>/02_geneformer_in_silico_screen/confounders/singlect_vs_jointpg_comparison.tsv
"""
import os
import re
import sys
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr
import gseapy as gp

# Section path configuration (GF_WORK, results folder, shipped data): see ../gf_paths.py.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gf_paths  # noqa: E402

FDR = 0.05
TBI = gf_paths.GF_WORK
JPG = os.path.join(gf_paths.MT_ISP, "garza_disease_binary_TBI_to_control",
                   "stats_delete", "by_cell_type")
JPGP = "isp_garza_disease_binary_TBI_to_control_delete"
MECHANO = gf_paths.MECHANO_GENES   # $GF_MECHANO/data/mechano_genes_list.txt
OUT_TSV = os.path.join(gf_paths.gf_results_dir("confounders"),
                       "singlect_vs_jointpg_comparison.tsv")

# cell-type name maps between the two pipelines
CTS = [("oligodendrocyte", "oligo", "oligodendrocyte"),
       ("astrocyte", "astro", "astrocyte"),
       ("microglia_macrophage", "microglia", "microglia_macrophage")]

PCG_TFS = {"EZH2", "EZH1", "SUZ12", "EED", "MTF2", "JARID2", "RING1B", "RNF2",
           "RING1", "BMI1", "CBX7", "CBX8", "PCGF2", "AEBP2"}


def apply_paper(df):
    return df[(df["goal_end_fdr"] < FDR) & (df["mean_shift"] > df["mean_shift"].mean())]


def fisher_set(elig, hit_ids, setids, idcol):
    bg = set(elig[idcol].dropna()); s = setids & bg
    a = len(s & hit_ids); b = len(s) - a; c = len(hit_ids - s); d = len(bg) - a - b - c
    if a + b == 0:
        return float("nan"), float("nan"), a
    orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return orr, p, a


# ---- membership sets ----
mech = set()
for ln in open(MECHANO):
    e = ln.split("|")[0].strip()
    if e.startswith("ENSG"):
        mech.add(e)
rea = gp.get_library("Reactome_2022", organism="Human")
ECM = set(rea[next(k for k in rea if "extracellular matrix organization" in k.lower())])

chea = gp.get_library("ChEA_2022", organism="Human")
enc = gp.get_library("ENCODE_TF_ChIP-seq_2015", organism="Human")
cnt = {}
for lib in (chea, enc):
    for term, genes in lib.items():
        if re.split(r"[ _]", term)[0].upper() in PCG_TFS:
            for g in set(genes):
                cnt[g] = cnt.get(g, 0) + 1
PCG = set(pd.Series(cnt)[pd.Series(cnt) >= 8].index)
print(f"sets: mechano {len(mech)} ENSG | ECM {len(ECM)} sym | consensus-PcG(K>=8) {len(PCG)} sym\n")

MECH_GENES = ["YAP1", "TEAD1", "WWTR1", "ACTN2", "CD44", "COL18A1", "TNC",
              "PTK2", "TLN2", "MRTFA", "DIAPH1", "CTNNA1", "AHNAK", "FN1"]

rows = []
for jpg_ct, tbi_ct, label in CTS:
    z = pd.read_csv(f"{TBI}/results/perturb_{tbi_ct}_paperQC_nGene1000/stats/"
                    f"{tbi_ct}_paperQC_nGene1000_filtered_high_confidence_gene_ranked.tsv", sep="\t")
    j = pd.read_csv(f"{JPG}/{JPGP}_{jpg_ct}_filtered_high_confidence_gene_ranked.tsv", sep="\t")
    hz = apply_paper(z); hj = apply_paper(j)
    zh = set(hz["ensembl_id"]); jh = set(hj["ensembl_id"])
    inter = zh & jh; uni = zh | jh
    jacc = len(inter) / len(uni) if uni else float("nan")
    # top-50 overlap by mean_shift
    z50 = set(hz.sort_values("mean_shift", ascending=False)["ensembl_id"].head(50))
    j50 = set(hj.sort_values("mean_shift", ascending=False)["ensembl_id"].head(50))
    top50 = len(z50 & j50)
    # Spearman on shared eligible genes (all eligible, not just hits)
    m = z[["ensembl_id", "mean_shift"]].merge(j[["ensembl_id", "mean_shift"]],
                                              on="ensembl_id", suffixes=("_z", "_j"))
    rho = spearmanr(m["mean_shift_z"], m["mean_shift_j"]).correlation if len(m) > 10 else float("nan")
    # enrichment in zero-shot hits (mirror of JointPG L1/L2)
    mo, mp, mk = fisher_set(z, zh, mech, "ensembl_id")
    eo, ep, ek = fisher_set(z, set(hz["gene_symbol"]), ECM, "gene_symbol")
    po, pp, pk = fisher_set(z, set(hz["gene_symbol"]), PCG, "gene_symbol")
    rows.append(dict(cell_type=label, zeroshot_hits=len(zh), jointpg_hits=len(jh),
                     shared=len(inter), jaccard=round(jacc, 3), top50_overlap=top50,
                     rho_eligible=round(rho, 3) if rho == rho else np.nan,
                     mech_OR=round(mo, 2), mech_p=mp, mech_k=mk,
                     ecm_OR=round(eo, 2), ecm_p=ep, ecm_k=ek,
                     pcg_OR=round(po, 2), pcg_p=pp, pcg_k=pk))
    # mechano-gene presence/rank in zero-shot
    print(f"=== {label}  (zero-shot {len(zh)} hits / JointPG {len(jh)} hits / shared {len(inter)}) ===")
    zr = z.set_index("gene_symbol")
    present = []
    for g in MECH_GENES:
        if g in zr.index:
            r = zr.loc[g]
            sig = "HIT" if (r["goal_end_fdr"] < FDR and r["mean_shift"] > z["mean_shift"].mean()) else "ns"
            present.append(f"{g}(rank{int(r['rank'])},fdr{r['goal_end_fdr']:.0e},{sig})")
    print("  mechano genes:", ", ".join(present) if present else "none in eligible set")

print("\n" + "=" * 110)
summ = pd.DataFrame(rows)
pd.set_option("display.width", 200, "display.max_columns", 40)
print(summ.to_string(index=False))
summ.to_csv(OUT_TSV, sep="\t", index=False)
print("\nwrote", OUT_TSV)
