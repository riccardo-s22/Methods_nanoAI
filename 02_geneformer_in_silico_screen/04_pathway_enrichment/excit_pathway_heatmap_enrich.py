#!/usr/bin/env python
"""Unified excitatory-neuron Enrichr across the 4 rescue arms used by the
pathway-level heatmap (MT_delete, MT_OE, FT_delete, FT_OE; all TBI->control).

All four arms run the SAME 20%-detection / top-100 / shared-background
methodology (mirrors mt_20pct_top100_enrichr.py) so the resulting per-arm
{label}_allterms.csv files are directly comparable cell-by-cell in the heatmap.

  - eligibility/background via the model-independent layer3 detection table
    (garza_qc4_de_by_celltype.tsv, excitatory_neuron det_TBI/det_ctrl),
  - foreground = FDR-significant hits (goal_end_fdr < 0.05, BH-corrected
    Wilcoxon; det>=20%; shift>0 toward goal) whose shift exceeds the AVERAGE
    shift OF THOSE SIGNIFICANT GENES (a sharper bar than the all-gene average,
    which is floored near zero by the tight shift distribution). No top-N cap.
  - background = ALL gene symbols detected >=20% in excitatory neurons,
  - libraries incl. MSigDB_Hallmark_2020.

Run with an arm label (excit only): MT_delete | FT_delete | FT_OE | MT_OE.

Locations (see ../gf_paths.py): the detection table and the MT ranked tables are
read below $GF_JOINTPG, the FT ranked tables below $GF_WORK (or, if absent there,
from data/02_geneformer_in_silico_screen/in_silico_perturbation/), and the output
goes to <results>/02_geneformer_in_silico_screen/pathway_heatmap/enrich/. The
detection table is not shipped, so this script cannot run from a fresh clone.
"""
import os
import sys
import pandas as pd
import gseapy as gp

# Section path configuration (GF_WORK, results folder, shipped data): see ../gf_paths.py.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gf_paths  # noqa: E402

JPG = gf_paths.GF_JOINTPG
TBI = gf_paths.GF_WORK
DET = gf_paths.MT_DETECTION_TABLE   # $GF_JOINTPG/results/isp/garza_disease_binary_TBI_to_control/mechano_analysis/layer3/garza_qc4_de_by_celltype.tsv
CT = "excitatory_neuron"
THRESH = 0.20
FDR = 0.05
LIBS = ["GO_Biological_Process_2026", "KEGG_2026",
        "Reactome_Pathways_2024", "MSigDB_Hallmark_2020"]

OUTROOT = os.path.join(gf_paths.GF_RESULTS, "pathway_heatmap", "enrich")

# label -> ranked-table path (all excitatory-neuron, all TBI->control rescue)
# MT: $GF_JOINTPG/results/isp/garza_disease_binary_<direction>/stats_<mode>/by_cell_type/
#     isp_garza_disease_binary_<direction>_<mode>_excitatory_neuron_unfiltered_all_gene_ranked.tsv
# FT: $GF_WORK/results/perturb_neurons_finetuned_<direction>/stats_<mode>/
#     excit_finetuned[_rescue]_<mode>_unfiltered_all_gene_ranked.tsv, or the same
#     table shipped in data/ when it is not below GF_WORK.
ARMS = {
    "MT_delete": gf_paths.mt_stats_table("TBI_to_control", "delete"),
    "MT_OE": gf_paths.mt_stats_table("TBI_to_control", "overexpress"),
    "FT_delete": gf_paths.ft_stats_table("rescue_delete"),
    "FT_OE": gf_paths.ft_stats_table("rescue_overexpress"),
    # MT off-diagonal direction (control->TBI = induction). goal=TBI, so
    # mean_shift>0 = toward TBI = a driver that INDUCES the TBI state.
    "MT_oe_c2t": gf_paths.mt_stats_table("control_to_TBI", "overexpress"),
    "MT_del_c2t": gf_paths.mt_stats_table("control_to_TBI", "delete"),
    # FT off-diagonal direction (control->TBI = induction; goal=TBI). Note the
    # induction tables drop the "rescue" infix in the filename.
    "FT_del_c2t": gf_paths.ft_stats_table("induction_delete"),
    "FT_oe_c2t": gf_paths.ft_stats_table("induction_overexpress"),
}


def clean(series):
    g = [str(x) for x in series.dropna().tolist()]
    return [x for x in g if x and x != "nan" and not x.startswith("ENSG")]


def run_arm(label):
    path = ARMS[label]
    if not os.path.exists(path):
        print(f"[{label}] ranked table MISSING -> {path}")
        return None
    unf = pd.read_csv(path, sep="\t")
    det = pd.read_csv(DET, sep="\t")
    det = det[det.cell_type == CT][["ensembl_id", "gene_symbol",
                                    "det_TBI", "det_ctrl"]].copy()
    det["det_max"] = det[["det_TBI", "det_ctrl"]].max(axis=1)
    back = clean(det[det["det_max"] >= THRESH]["gene_symbol"])
    m = unf.merge(det[["ensembl_id", "det_TBI", "det_ctrl", "det_max"]],
                  on="ensembl_id", how="inner")
    # FDR-significant, det>=20%, toward goal (shift>0)
    sig = m[(m["det_max"] >= THRESH) & (m["mean_shift"] > 0)
            & (m["goal_end_fdr"] < FDR)].copy()
    # bar = average shift OF the significant genes (sharper than all-gene avg)
    avg_sig = float(sig["mean_shift"].mean()) if len(sig) else 0.0
    fg_df = sig[sig["mean_shift"] > avg_sig].sort_values("mean_shift", ascending=False)
    fg = clean(fg_df["gene_symbol"])
    print(f"\n##### {label} #####")
    print(f"  perturbed&matched: {len(m):,} | FDR-sig(det>=20%,shift>0): "
          f"{len(sig):,} | avg_shift_of_sig={avg_sig:.4g} | "
          f"HITS (sig & shift>sig-avg): {len(fg_df):,} | background det>=20%: {len(back):,}")
    print(f"  top15 by shift: {', '.join(fg[:15])}")
    outdir = f"{OUTROOT}/{label}"
    os.makedirs(outdir, exist_ok=True)
    keep = [c for c in ["rank", "gene_symbol", "ensembl_id", "mean_shift",
                        "goal_end_fdr", "n_cells_perturbed", "det_TBI",
                        "det_ctrl"] if c in fg_df.columns]
    fg_df[keep].to_csv(f"{outdir}/{label}_hits_fdr05.tsv",
                       sep="\t", index=False)
    if len(fg) < 20:
        print(f"  SKIP enrichment: only {len(fg)} eligible genes")
        return None
    enr = gp.enrichr(gene_list=fg, gene_sets=LIBS, background=back,
                     outdir=outdir, no_plot=True, verbose=False)
    res = enr.results
    res.to_csv(f"{outdir}/{label}_enrichr_allterms.csv", index=False)
    sig = res[res["Adjusted P-value"] < 0.05]
    print(f"  Enrichr sig(adjP<0.05): {len(sig)} terms  "
          f"(min adjP={res['Adjusted P-value'].min():.1e})")
    for lib in LIBS:
        nsg = int((sig["Gene_set"] == lib).sum())
        print(f"    {lib}: {nsg} sig")
    return res


if __name__ == "__main__":
    if not os.path.exists(DET):
        raise SystemExit(
            f"ERROR: detection table not found: {DET}\n"
            "It defines the 7,374-gene eligibility universe and is not shipped in "
            "data/. Set GF_JOINTPG to a tree that holds it (see ../gf_paths.py).")
    os.makedirs(OUTROOT, exist_ok=True)
    todo = sys.argv[1:] if len(sys.argv) > 1 else list(ARMS)
    for label in todo:
        if label not in ARMS:
            print(f"unknown arm '{label}'; choices: {list(ARMS)}")
            continue
        run_arm(label)
    print(f"\nDONE -> {OUTROOT}/<label>/<label>_enrichr_allterms.csv")
