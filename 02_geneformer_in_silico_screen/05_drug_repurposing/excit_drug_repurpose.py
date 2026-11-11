#!/usr/bin/env python
"""Drug-repurposing screen: compounds that may PREVENT excitatory neurons from
shifting toward the TBI state, built on the ISP driver gene lists.

Direction-aware union (per the connectivity-map logic):
  INHIBIT targets  (drug should LOWER the gene; query drug DOWN signatures):
      induction-OE hits   (overexpress -> TBI ; raising is pathogenic)
    U rescue-delete hits  (delete      -> ctrl; lowering is therapeutic)
  ACTIVATE targets (drug should RAISE the gene; query drug UP signatures):
      induction-delete hits (delete     -> TBI ; the gene is protective)
    U rescue-OE hits        (overexpress-> ctrl; raising is protective)

A drug that DOWN-regulates INHIBIT targets AND UP-regulates ACTIVATE targets
opposes the TBI shift on both axes -> top candidate. Genes that land in BOTH
target lists are direction-ambiguous and dropped. Genes in the bidirectional
INTERSECTION (induction AND rescue agree) are the highest-confidence targets and
are reported separately.

Model = FT (excit-specific fine-tuned). Background = genes det>=20% in excit.

Locations (see ../gf_paths.py): the detection table and the MT ranked tables are
read below $GF_JOINTPG, the FT ranked tables below $GF_WORK (or, if absent there,
from data/02_geneformer_in_silico_screen/in_silico_perturbation/), and the output
goes to <results>/02_geneformer_in_silico_screen/drug_repurpose/<model>/. The
detection table is not shipped, so this script cannot run from a fresh clone.
"""
import os
import sys
import pandas as pd
import gseapy as gp

# Section path configuration (GF_WORK, results folder, shipped data): see ../gf_paths.py.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gf_paths  # noqa: E402

TBI = gf_paths.GF_WORK
JPG = gf_paths.GF_JOINTPG
ENR = gf_paths.analysis_input("pathway_heatmap", "enrich")
DET = gf_paths.MT_DETECTION_TABLE   # $GF_JOINTPG/results/isp/garza_disease_binary_TBI_to_control/mechano_analysis/layer3/garza_qc4_de_by_celltype.tsv
CT, THRESH, FDR = "excitatory_neuron", 0.20, 0.05
N_CELLS_MIN = 50      # min perturbed cells per gene
SHIFT_PCT = 0.75      # keep shift > this per-arm quantile of the significant set
OUT = os.path.join(gf_paths.GF_RESULTS, "drug_repurpose")

# raw ranked-table paths per arm label (same sources as excit_pathway_heatmap_enrich.py)
# MT: $GF_JOINTPG/results/isp/garza_disease_binary_<direction>/stats_<mode>/by_cell_type/
#     isp_garza_disease_binary_<direction>_<mode>_excitatory_neuron_unfiltered_all_gene_ranked.tsv
# FT: $GF_WORK/results/perturb_neurons_finetuned_<direction>/stats_<mode>/
#     excit_finetuned[_rescue]_<mode>_unfiltered_all_gene_ranked.tsv, or the same
#     table shipped in data/ when it is not below GF_WORK.
RANKED = {
    "MT_delete": gf_paths.mt_stats_table("TBI_to_control", "delete"),
    "MT_OE": gf_paths.mt_stats_table("TBI_to_control", "overexpress"),
    "MT_oe_c2t": gf_paths.mt_stats_table("control_to_TBI", "overexpress"),
    "MT_del_c2t": gf_paths.mt_stats_table("control_to_TBI", "delete"),
    "FT_delete": gf_paths.ft_stats_table("rescue_delete"),
    "FT_OE": gf_paths.ft_stats_table("rescue_overexpress"),
    "FT_del_c2t": gf_paths.ft_stats_table("induction_delete"),
    "FT_oe_c2t": gf_paths.ft_stats_table("induction_overexpress"),
}
_DET_CACHE = None

MODEL = sys.argv[1] if len(sys.argv) > 1 else "FT"   # FT | MT
MODE = sys.argv[2] if len(sys.argv) > 2 else "union"  # union | core

# LINCS/GEO library aliases -> canonical drug name (dedupe naming variants)
ALIAS = {"plx4032": "vemurafenib", "plx-4032": "vemurafenib",
         " vorinostat": "vorinostat"}
# non-therapeutic perturbagens / industrial toxicants / chemo to flag-out of the
# headline (library noise, not repurposing candidates)
ARTIFACTS = {"hypochlorous", "phosgene", "monuron", "laccaic", "mt19c",
             "2,3,7,8-tetrachlorodibenzo-p-dioxin", "lipopolysaccharide",
             "halothane", "cisplatin", "carboplatin", "oxaliplatin",
             "paraquat", "rotenone", "tunicamycin", "thapsigargin",
             "hydrogen", "ozone", "formaldehyde", "benzo", "acrolein"}
A = {  # arm -> hits file label, by model
    "ind_oe":  f"{MODEL}_oe_c2t",    # overexpress -> TBI  -> INHIBIT
    "res_del": f"{MODEL}_delete",    # delete      -> ctrl -> INHIBIT
    "ind_del": f"{MODEL}_del_c2t",   # delete      -> TBI  -> ACTIVATE
    "res_oe":  f"{MODEL}_OE",        # overexpress -> ctrl -> ACTIVATE
}
DOWN_LIBS = ["LINCS_L1000_Chem_Pert_down", "Drug_Perturbations_from_GEO_down"]
UP_LIBS = ["LINCS_L1000_Chem_Pert_up", "Drug_Perturbations_from_GEO_up"]


def _det():
    global _DET_CACHE
    if _DET_CACHE is None:
        d = pd.read_csv(DET, sep="\t")
        d = d[d.cell_type == CT].copy()
        d["det_max"] = d[["det_TBI", "det_ctrl"]].max(axis=1)
        _DET_CACHE = d
    return _DET_CACHE


def _clean(series):
    g = series.dropna().astype(str)
    return set(x for x in g if x and x != "nan" and not x.startswith("ENSG"))


def hits(label):
    """Foreground = det>=20%, shift>0 toward goal, FDR-sig, n_cells>=50, AND
    shift above the per-arm 75th percentile of the significant set."""
    p = RANKED[label]
    if not os.path.exists(p):
        print(f"  [{label}] ranked table MISSING -> empty set ({p})")
        return set()
    unf = pd.read_csv(p, sep="\t")
    m = unf.merge(_det()[["ensembl_id", "det_max"]], on="ensembl_id", how="inner")
    sig = m[(m["det_max"] >= THRESH) & (m["mean_shift"] > 0)
            & (m["goal_end_fdr"] < FDR) & (m["n_cells_perturbed"] >= N_CELLS_MIN)].copy()
    if len(sig) == 0:
        return set()
    bar = sig["mean_shift"].quantile(SHIFT_PCT)
    fg = sig[sig["mean_shift"] > bar]
    return _clean(fg["gene_symbol"])


def background():
    return sorted(_clean(_det()[_det()["det_max"] >= THRESH]["gene_symbol"]))


def enrich(genes, libs, tag, back):
    genes = sorted(genes)
    if len(genes) < 15:
        print(f"  [{tag}] only {len(genes)} genes -> skip enrichment")
        return pd.DataFrame()
    outdir = f"{OUT}/{MODEL}/{tag}"
    os.makedirs(outdir, exist_ok=True)
    enr = gp.enrichr(gene_list=genes, gene_sets=libs, background=back,
                     outdir=outdir, no_plot=True, verbose=False)
    r = enr.results.copy()
    r.to_csv(f"{outdir}/{tag}_allterms.csv", index=False)
    return r


def drug_name(term):
    # LINCS terms look like "withaferin-a MCF7 ..."; GEO like "drugname ...".
    d = str(term).split(" ")[0].lower().strip()
    return ALIAS.get(d, d)


if __name__ == "__main__":
    if not os.path.exists(DET):
        raise SystemExit(
            f"ERROR: detection table not found: {DET}\n"
            "It defines the 7,374-gene eligibility universe and is not shipped in "
            "data/. Set GF_JOINTPG to a tree that holds it (see ../gf_paths.py).")
    os.makedirs(f"{OUT}/{MODEL}", exist_ok=True)
    print(f"===== model {MODEL} =====")
    s = {k: hits(v) for k, v in A.items()}
    for k, v in A.items():
        print(f"  {k:8s} ({v}): {len(s[k])} genes")

    inhibit = s["ind_oe"] | s["res_del"]
    activate = s["ind_del"] | s["res_oe"]
    conflict = inhibit & activate
    inhibit -= conflict
    activate -= conflict
    inh_core = (s["ind_oe"] & s["res_del"]) - conflict   # bidirectional concordant
    act_core = (s["ind_del"] & s["res_oe"]) - conflict
    print(f"\n  INHIBIT targets (drug down): {len(inhibit)}  "
          f"[bidirectional core: {len(inh_core)}]")
    print(f"  ACTIVATE targets (drug up):  {len(activate)}  "
          f"[bidirectional core: {len(act_core)}]")
    print(f"  dropped direction-conflict genes: {len(conflict)}")
    print(f"  INHIBIT core genes: {', '.join(sorted(inh_core)[:30])}")
    print(f"  ACTIVATE core genes: {', '.join(sorted(act_core)[:30])}")

    if MODE == "core":
        inhibit, activate = inh_core, act_core
        print(f"\n  MODE=core -> using bidirectional-concordant targets only "
              f"(INHIBIT {len(inhibit)}, ACTIVATE {len(activate)})")
    else:
        print(f"\n  MODE=union -> using full direction-aware union")

    back = background()
    print(f"  background det>=20%: {len(back)}")

    print("\n--- enrichment: INHIBIT targets vs drug DOWN signatures ---")
    dn = enrich(inhibit, DOWN_LIBS, f"inhibit_down_{MODE}", back)
    print("--- enrichment: ACTIVATE targets vs drug UP signatures ---")
    up = enrich(activate, UP_LIBS, f"activate_up_{MODE}", back)

    # combine: drug -> summed -log10(adjP) across both directions
    import numpy as np
    rows = []
    for r, direction in [(dn, "down/inhibit"), (up, "up/activate")]:
        if r.empty:
            continue
        rr = r[r["Adjusted P-value"] < 0.05].copy()
        rr["drug"] = rr["Term"].map(drug_name)
        rr["mlog10adjp"] = -np.log10(rr["Adjusted P-value"].clip(lower=1e-300))
        rr["direction"] = direction
        rows.append(rr[["drug", "direction", "Term", "Gene_set",
                        "Adjusted P-value", "Odds Ratio", "mlog10adjp", "Genes"]])
    if rows:
        allr = pd.concat(rows, ignore_index=True)
        allr.to_csv(f"{OUT}/{MODEL}/all_sig_drug_terms.tsv", sep="\t", index=False)
        # per-drug: which directions it hits, summed signal, n terms
        agg = (allr.groupby("drug")
               .agg(n_terms=("Term", "size"),
                    directions=("direction", lambda x: "+".join(sorted(set(x)))),
                    sum_signal=("mlog10adjp", "sum"),
                    best_adjp=("Adjusted P-value", "min"),
                    max_OR=("Odds Ratio", "max"))
               .reset_index())
        # bonus weight for drugs hitting BOTH axes (down+up)
        agg["both_axes"] = agg["directions"].str.contains(r"\+")
        agg["artifact"] = agg["drug"].isin(ARTIFACTS)
        agg = agg.sort_values(["artifact", "both_axes", "sum_signal"],
                              ascending=[True, False, False]).reset_index(drop=True)
        out = f"{OUT}/{MODEL}/drug_candidates_ranked_{MODE}.tsv"
        agg.to_csv(out, sep="\t", index=False)
        clean = agg[~agg["artifact"]]
        print(f"\n  sig drug terms: {len(allr)} | unique drugs: {len(agg)} | "
              f"both-axes: {int(agg['both_axes'].sum())} | "
              f"artifacts flagged: {int(agg['artifact'].sum())}")
        print("\n  TOP 25 candidate drugs (artifacts excluded):")
        with pd.option_context("display.max_rows", 30, "display.width", 160):
            print(clean.head(25).to_string(index=False))
        print(f"\nwrote {out}")
    else:
        print("\n  no significant drug terms in either direction")
