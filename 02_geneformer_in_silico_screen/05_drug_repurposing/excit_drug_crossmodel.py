#!/usr/bin/env python
"""Cross-model drug consensus: keep compounds that oppose the excit TBI shift in
BOTH the FT (excit-specific) and MT models. Reads the per-model union ranked
tables, drops flagged library artifacts, intersects on drug, and marks drugs
that also survive the stringent FT bidirectional-concordant-core list.
"""
import pandas as pd

OUT = "/work/pi_cfallini_uri_edu/Geneformer_TBI/results/drug_repurpose"


def load(model, mode):
    df = pd.read_csv(f"{OUT}/{model}/drug_candidates_ranked_{mode}.tsv", sep="\t")
    return df[~df["artifact"]].copy()


ft = load("FT", "union")
mt = load("MT", "union")
ft_core = set(load("FT", "core")["drug"])

m = ft.merge(mt, on="drug", suffixes=("_FT", "_MT"))
m["combined_signal"] = m["sum_signal_FT"] + m["sum_signal_MT"]
m["in_FT_core"] = m["drug"].isin(ft_core)
m["both_axes_either"] = m["both_axes_FT"] | m["both_axes_MT"]
m = m.sort_values(["in_FT_core", "both_axes_either", "combined_signal"],
                  ascending=[False, False, False]).reset_index(drop=True)

cols = ["drug", "in_FT_core", "both_axes_either", "combined_signal",
        "directions_FT", "directions_MT", "best_adjp_FT", "best_adjp_MT",
        "n_terms_FT", "n_terms_MT"]
m[cols].to_csv(f"{OUT}/crossmodel_consensus.tsv", sep="\t", index=False)
print(f"FT union drugs (clean): {len(ft)} | MT union drugs (clean): {len(mt)}")
print(f"consensus (both models): {len(m)}  | also in FT-core: {int(m['in_FT_core'].sum())}")
with pd.option_context("display.width", 170, "display.max_rows", 60):
    print(m[cols].to_string(index=False))
print(f"\nwrote {OUT}/crossmodel_consensus.tsv")
