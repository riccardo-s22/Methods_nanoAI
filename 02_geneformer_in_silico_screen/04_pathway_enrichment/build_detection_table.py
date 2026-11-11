#!/usr/bin/env python
"""Rebuild the 20%-detection table read by the pathway and drug scripts.

NOT part of the original run. In the published analysis the table
garza_qc4_de_by_celltype.tsv was written by the MT project, whose code and
files were lost; it is not shipped in data/ and no original script writes it.
excit_pathway_heatmap_enrich.py, ../05_drug_repurposing/excit_drug_repurpose.py
and ../figures/heatmap_excit_finetuned_cosine_shifts.py cannot run without it.
This helper was added so that they can run once 01_dataset_preparation has
produced the excitatory-neuron AnnData. It selects no genes and runs no test.

What it computes is the eligibility definition of the Methods: for every gene,
the fraction of excitatory nuclei of each disease state with a nonzero raw
count. The expression is the one confound_check_sex_age_drivers.py already uses
for its detection universe.

Only what the downstream scripts read is written: the excitatory_neuron rows
and the columns cell_type, ensembl_id, gene_symbol, det_TBI, det_ctrl. The
original table also held other cell types and other columns, which nothing in
this repository reads.

Checked against the original run: the per-arm driver lists of that run carry
det_TBI and det_ctrl for their genes (6,093 distinct genes over the eight FT
and MT arms); the rebuilt values agree with all of them to within 1.2e-16 (one
unit in the last place), no gene lies within 1e-9 of the 0.20 threshold, the
table gives 7,374 genes detected in >=20% of either group, and
excit_pathway_heatmap_enrich.py run on it returns driver lists byte-identical
to the original ones for the four FT arms. Genes outside those lists could not
be compared, because the original table is not available.

Input   $GF_WORK/h5ad/celltype_subsets_paperQC_nGene1000/excitatory_neurons_raw_geneformer_ready.h5ad
Output  $GF_JOINTPG/results/isp/garza_disease_binary_TBI_to_control/mechano_analysis/layer3/garza_qc4_de_by_celltype.tsv
        (gf_paths.MT_DETECTION_TABLE, where the downstream scripts look for it).
        An existing file is never overwritten.

Usage:  python 04_pathway_enrichment/build_detection_table.py     (needs anndata)
"""
import os
import sys

import anndata
import numpy as np
import pandas as pd

# Section path configuration (GF_WORK, GF_JOINTPG): see ../gf_paths.py.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gf_paths  # noqa: E402

H5AD = os.path.join(gf_paths.GF_WORK, "h5ad", "celltype_subsets_paperQC_nGene1000",
                    "excitatory_neurons_raw_geneformer_ready.h5ad")
OUT = gf_paths.MT_DETECTION_TABLE
CT = "excitatory_neuron"
THRESH = 0.20   # reported only; the downstream scripts apply the threshold themselves


def main():
    if os.path.exists(OUT):
        raise SystemExit(f"detection table already exists, not overwritten: {OUT}")
    if not os.path.exists(H5AD):
        raise SystemExit(f"ERROR: excitatory-neuron AnnData not found: {H5AD}\n"
                         "It is written by 01_dataset_preparation (see the section README).")
    a = anndata.read_h5ad(H5AD)
    X = a.X.tocsr()
    det = {c: np.asarray((X[(a.obs.condition == c).values] > 0).mean(0)).ravel()
           for c in ["TBI", "control"]}
    table = pd.DataFrame({
        "cell_type": CT,
        "ensembl_id": a.var["ensembl_id"].astype(str).values,
        "gene_symbol": a.var["gene_symbol"].astype(str).values,
        "det_TBI": det["TBI"],
        "det_ctrl": det["control"],
    })
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    table.to_csv(OUT, sep="\t", index=False)
    n = a.obs.condition.value_counts()
    eligible = table[table[["det_TBI", "det_ctrl"]].max(axis=1) >= THRESH]
    print(f"nuclei: {a.n_obs:,} ({int(n.get('TBI', 0)):,} TBI, {int(n.get('control', 0)):,} control); "
          f"genes: {a.n_vars:,}")
    print(f"genes detected in >={THRESH:.0%} of either group: {len(eligible):,} "
          f"({eligible.gene_symbol.nunique():,} distinct symbols)")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
