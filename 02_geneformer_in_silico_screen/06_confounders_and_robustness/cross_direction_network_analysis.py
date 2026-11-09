#!/usr/bin/env python
"""Cross-direction network relationship analysis for fine-tuned excitatory
neuron in-silico perturbation screens.

v2 — four datasets (all on the same 2-class fine-tuned TBI/control axis,
excitatory neurons):
  D1 = INDUCER-delete       (control -> TBI by deletion)
  D2 = INDUCER-overexpress  (control -> TBI by overexpression)
  D3 = RESCUE-delete        (TBI -> control by deletion)
  D4 = RESCUE-overexpress   (TBI -> control by overexpression)   [new in v2]

Three views of the relationship:
  1. Gene-list overlap with hypergeometric significance vs. tested universe
     (mode symmetry, mode complementarity, dose-response candidates)
  2. Pathway-level convergence across Enrichr libraries
  3. STRING PPI subgraph of the union top-N, with shared-neighbor counts
     and connected-component clustering

Outputs land under
  /work/pi_cfallini_uri_edu/Geneformer_TBI/results/cross_direction_network_v2/
"""
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import requests
from scipy.stats import hypergeom

BASE = Path("/work/pi_cfallini_uri_edu/Geneformer_TBI/results")
OUT  = BASE / "cross_direction_network_v2"
OUT.mkdir(parents=True, exist_ok=True)

DATASETS = {
    "INDUCER_del": {
        "label": "control -> TBI  (delete)",
        "bio":   BASE / "perturb_neurons_finetuned_control_to_tbi/stats_delete/excit_finetuned_delete_filtered_biological_candidate_ranked.tsv",
        "all":   BASE / "perturb_neurons_finetuned_control_to_tbi/stats_delete/excit_finetuned_delete_unfiltered_all_gene_ranked.tsv",
        "enr":   BASE / "perturb_neurons_finetuned_control_to_tbi/enrichment_delete",
    },
    "INDUCER_ovr": {
        "label": "control -> TBI  (overexpress)",
        "bio":   BASE / "perturb_neurons_finetuned_control_to_tbi/stats_overexpress/excit_finetuned_overexpress_filtered_biological_candidate_ranked.tsv",
        "all":   BASE / "perturb_neurons_finetuned_control_to_tbi/stats_overexpress/excit_finetuned_overexpress_unfiltered_all_gene_ranked.tsv",
        "enr":   BASE / "perturb_neurons_finetuned_control_to_tbi/enrichment_overexpress",
    },
    "RESCUE_del": {
        "label": "TBI -> control  (delete)",
        "bio":   BASE / "perturb_neurons_finetuned_tbi_to_control/stats_delete/excit_finetuned_rescue_delete_filtered_biological_candidate_ranked.tsv",
        "all":   BASE / "perturb_neurons_finetuned_tbi_to_control/stats_delete/excit_finetuned_rescue_delete_unfiltered_all_gene_ranked.tsv",
        "enr":   BASE / "perturb_neurons_finetuned_tbi_to_control/enrichment_delete",
    },
    "RESCUE_ovr": {
        "label": "TBI -> control  (overexpress)",
        "bio":   BASE / "perturb_neurons_finetuned_tbi_to_control/stats_overexpress/excit_finetuned_rescue_overexpress_filtered_biological_candidate_ranked.tsv",
        "all":   BASE / "perturb_neurons_finetuned_tbi_to_control/stats_overexpress/excit_finetuned_rescue_overexpress_unfiltered_all_gene_ranked.tsv",
        "enr":   BASE / "perturb_neurons_finetuned_tbi_to_control/enrichment_overexpress",
    },
}

TOP_N        = 200      # how many top hits per dataset for overlap + STRING
STRING_SCORE = 700      # high-confidence STRING interaction cutoff (0-1000)
ENRICHMENT_FDR = 0.05


# ============================================================
# 1. Load datasets
# ============================================================

print("=" * 60)
print("Cross-direction network relationship analysis")
print(f"TOP_N per dataset: {TOP_N}    STRING score cutoff: {STRING_SCORE}")
print("=" * 60)

data = {}
for key, cfg in DATASETS.items():
    bio = pd.read_csv(cfg["bio"], sep="\t")
    uni = pd.read_csv(cfg["all"], sep="\t", usecols=["gene_symbol"])
    bio = bio[bio["sig"] == 1].copy().reset_index(drop=True)
    bio_sorted = bio.sort_values("mean_shift", ascending=False)
    top = bio_sorted.head(TOP_N)["gene_symbol"].astype(str).tolist()
    data[key] = {
        "label":     cfg["label"],
        "top":       set(top),
        "top_list":  top,
        "bio_genes": set(bio["gene_symbol"].astype(str)),
        "universe":  set(uni["gene_symbol"].astype(str)),
        "bio_df":    bio,
        "enr_dir":   cfg["enr"],
    }
    print(f"  {key:12s}  ({cfg['label']})")
    print(f"      universe: {len(data[key]['universe']):>6,}  bio_sig: {len(data[key]['bio_genes']):>6,}  top_{TOP_N}: {len(top)}")


# ============================================================
# 2. Pairwise overlap + hypergeometric
# ============================================================
print("\n" + "=" * 60)
print("2. Pairwise top-{} overlap with hypergeometric significance".format(TOP_N))
print("=" * 60)

pair_rows = []
intersections = {}
keys = list(data.keys())
for i, a in enumerate(keys):
    for b in keys[i + 1:]:
        A, B = data[a]["top"], data[b]["top"]
        U    = data[a]["universe"] | data[b]["universe"]
        inter = A & B
        union = A | B
        jacc  = len(inter) / len(union) if union else 0.0
        # Hypergeometric: prob of drawing |A & B| or more shared between two
        # random TOP_N-size samples from U.
        M = len(U)
        n = len(A)
        N = len(B)
        k = len(inter)
        # P(X >= k) where X ~ Hypergeometric(M, n, N)
        pval = hypergeom.sf(k - 1, M, n, N) if k > 0 else 1.0
        fold = (k / N) / (n / M) if (N > 0 and n > 0 and M > 0) else 0.0
        pair_rows.append(dict(
            set_A=a, set_B=b, A_size=len(A), B_size=len(B),
            universe=M, intersection=k, jaccard=round(jacc, 4),
            fold_enrichment=round(fold, 2), hypergeometric_p=pval
        ))
        intersections[(a, b)] = sorted(inter)
        print(f"  {a:12s} ∩ {b:12s}  k={k:>3}  jacc={jacc:.3f}  fold={fold:.2f}x  p={pval:.2e}")

pd.DataFrame(pair_rows).to_csv(OUT / "01_pairwise_overlap.tsv", sep="\t", index=False)

# Quadruple intersection (genes in top-N of every direction — strongest signal)
quad = data["INDUCER_del"]["top"] & data["INDUCER_ovr"]["top"] & data["RESCUE_del"]["top"] & data["RESCUE_ovr"]["top"]
print(f"\n  Quadruple intersection (all 4 top-{TOP_N}): {len(quad)} gene(s)")
if quad:
    print(f"    {sorted(quad)}")

# Each 3-of-4 triple — meaningful biologically since each leaves out one direction
triples = {
    "INDUCER_del ∩ INDUCER_ovr ∩ RESCUE_del": data["INDUCER_del"]["top"] & data["INDUCER_ovr"]["top"] & data["RESCUE_del"]["top"],
    "INDUCER_del ∩ INDUCER_ovr ∩ RESCUE_ovr": data["INDUCER_del"]["top"] & data["INDUCER_ovr"]["top"] & data["RESCUE_ovr"]["top"],
    "INDUCER_del ∩ RESCUE_del ∩ RESCUE_ovr":  data["INDUCER_del"]["top"] & data["RESCUE_del"]["top"]  & data["RESCUE_ovr"]["top"],
    "INDUCER_ovr ∩ RESCUE_del ∩ RESCUE_ovr":  data["INDUCER_ovr"]["top"] & data["RESCUE_del"]["top"]  & data["RESCUE_ovr"]["top"],
}
print(f"\n  3-of-4 intersections:")
for name, s in triples.items():
    print(f"    {name:45s} {len(s):>3} genes")

# Save intersection gene lists with context
inter_rows = []
for (a, b), genes in intersections.items():
    for g in genes:
        rec = {"set_A": a, "set_B": b, "gene": g}
        for ds in (a, b):
            row = data[ds]["bio_df"]
            hit = row[row["gene_symbol"] == g]
            if not hit.empty:
                r = hit.iloc[0]
                rec[f"{ds}_rank"] = int(r["rank"])
                rec[f"{ds}_shift"] = round(float(r["mean_shift"]), 5)
                rec[f"{ds}_n_donors"] = int(r["n_donors_contributing"])
                rec[f"{ds}_consistency"] = round(float(r["donor_consistency_score"]), 3)
        inter_rows.append(rec)
pd.DataFrame(inter_rows).to_csv(OUT / "02_intersection_gene_table.tsv", sep="\t", index=False)


# ============================================================
# 3. Directional categorization
# ============================================================
print("\n" + "=" * 60)
print("3. Reversibility / dose-response categorization")
print("=" * 60)

# Each set is "top-N genes whose perturbation pushes the start state toward the goal state".
# Useful categories (all six pairwise overlaps are interpretable):
#   REVERSIBILITY axes — clean drivers; perturbation in one direction induces TBI,
#   the inverse perturbation rescues TBI:
#     (a) AXIS 1:  INDUCER_ovr ∩ RESCUE_del   (X high drives TBI; deleting X rescues)
#     (b) AXIS 2:  INDUCER_del ∩ RESCUE_ovr   (X low drives TBI; overexpressing X rescues)
#   MODE-SYMMETRY axes — same mode pushes both directions; X has state-dependent role:
#     (c) DEL_SYM: INDUCER_del ∩ RESCUE_del   (deleting X converges both starting states)
#     (d) OVR_SYM: INDUCER_ovr ∩ RESCUE_ovr   (overexpressing X converges both)
#   DOSE-RESPONSE within direction — both delete and overexpress push toward the same goal,
#   meaning that direction's start state is sensitive to either-sided dysregulation of X:
#     (e) DOSE_INDUCER: INDUCER_del ∩ INDUCER_ovr
#     (f) DOSE_RESCUE:  RESCUE_del  ∩ RESCUE_ovr

cats = {
    "reversibility AXIS 1 (INDUCER_ovr ∩ RESCUE_del)": data["INDUCER_ovr"]["top"] & data["RESCUE_del"]["top"],
    "reversibility AXIS 2 (INDUCER_del ∩ RESCUE_ovr)": data["INDUCER_del"]["top"] & data["RESCUE_ovr"]["top"],
    "del-mode symmetry   (INDUCER_del ∩ RESCUE_del)":  data["INDUCER_del"]["top"] & data["RESCUE_del"]["top"],
    "ovr-mode symmetry   (INDUCER_ovr ∩ RESCUE_ovr)":  data["INDUCER_ovr"]["top"] & data["RESCUE_ovr"]["top"],
    "dose-response INDUCER (INDUCER_del ∩ INDUCER_ovr)": data["INDUCER_del"]["top"] & data["INDUCER_ovr"]["top"],
    "dose-response RESCUE  (RESCUE_del  ∩ RESCUE_ovr)":  data["RESCUE_del"]["top"]  & data["RESCUE_ovr"]["top"],
}
for name, s in cats.items():
    print(f"  {name:50s} {len(s):>3} genes")

# ============================================================
# 4. Pathway-level convergence
# ============================================================
print("\n" + "=" * 60)
print(f"4. Pathway-level convergence (FDR < {ENRICHMENT_FDR})")
print("=" * 60)

LIBS = ["GO_Biological_Process_2023", "Reactome_2022", "KEGG_2021_Human", "MSigDB_Hallmark_2020"]
TIERS = ["high_tier"]  # focus on high-tier (top biological priority) for cross-comparison

pathway_hits = defaultdict(dict)  # term -> {ds: (fdr, overlap, genes)}
for ds_key, ds in data.items():
    enr_dir = ds["enr_dir"]
    if not enr_dir.exists():
        print(f"  [warn] no enrichment dir for {ds_key}: {enr_dir}")
        continue
    for tier in TIERS:
        for lib in LIBS:
            f = enr_dir / tier / f"{lib}.tsv"
            if not f.exists():
                continue
            df = pd.read_csv(f, sep="\t")
            fdr_col = next((c for c in df.columns if "Adjusted" in c or "adjusted" in c), None)
            term_col = "Term" if "Term" in df.columns else df.columns[0]
            ovl_col = "Overlap" if "Overlap" in df.columns else None
            gene_col = "Genes" if "Genes" in df.columns else None
            if fdr_col is None:
                continue
            sig = df[df[fdr_col] < ENRICHMENT_FDR]
            for _, row in sig.iterrows():
                key = (lib, row[term_col])
                pathway_hits[key][ds_key] = (
                    float(row[fdr_col]),
                    row.get(ovl_col, ""),
                    row.get(gene_col, ""),
                )

# Find pathways significant in >= 2 datasets
convergent = []
for (lib, term), ds_map in pathway_hits.items():
    if len(ds_map) >= 2:
        rec = {"library": lib, "term": term, "n_datasets": len(ds_map)}
        for k in DATASETS:
            if k in ds_map:
                fdr, ovl, genes = ds_map[k]
                rec[f"{k}_fdr"] = fdr
                rec[f"{k}_overlap"] = ovl
                rec[f"{k}_genes"] = genes
            else:
                rec[f"{k}_fdr"] = None
        convergent.append(rec)

conv_df = pd.DataFrame(convergent).sort_values("n_datasets", ascending=False) if convergent else pd.DataFrame()
conv_df.to_csv(OUT / "03_pathway_convergence.tsv", sep="\t", index=False)
print(f"  Pathways sig in >= 2 datasets: {len(conv_df)}")
if not conv_df.empty:
    print(f"  Pathways sig in all 4 datasets: {(conv_df['n_datasets'] == 4).sum()}")
    print(f"  Pathways sig in >= 3 datasets:  {(conv_df['n_datasets'] >= 3).sum()}")
    print("\n  Top 12 convergent pathways:")
    show_cols = ["n_datasets", "library", "term"]
    print(conv_df[show_cols].head(12).to_string(index=False))


# ============================================================
# 5. STRING PPI subgraph
# ============================================================
print("\n" + "=" * 60)
print(f"5. STRING PPI (high-confidence ≥{STRING_SCORE}) on union top-{TOP_N}")
print("=" * 60)

union_top = sorted(set().union(*(data[k]["top"] for k in DATASETS)))
print(f"  Union top-{TOP_N} across all {len(DATASETS)} datasets: {len(union_top)} unique genes")

# STRING API — POST identifiers in chunks (STRING limit is 2000 ids per call).
STRING_URL = "https://string-db.org/api/json/network"
session = requests.Session()
session.headers.update({"caller_identity": "geneformer_tbi_riccardo-s_uri_edu"})

print(f"  Querying STRING ({STRING_URL}) ...")
chunk = 1500
edges = []
for i in range(0, len(union_top), chunk):
    batch = union_top[i:i + chunk]
    params = {
        "identifiers": "\r".join(batch),
        "species":     9606,
        "required_score": STRING_SCORE,
        "network_type":   "functional",
    }
    try:
        r = session.post(STRING_URL, data=params, timeout=90)
        r.raise_for_status()
        batch_edges = r.json()
        edges.extend(batch_edges)
        print(f"    batch {i // chunk + 1}: {len(batch)} ids -> {len(batch_edges)} edges")
        time.sleep(1.0)
    except Exception as e:
        print(f"    [error] STRING query failed: {type(e).__name__}: {e}")
        sys.exit(2)

# Dedupe edges (STRING returns both directions for an undirected edge)
edge_pairs = set()
deduped = []
for e in edges:
    a, b = e["preferredName_A"], e["preferredName_B"]
    key = tuple(sorted([a, b]))
    if key in edge_pairs:
        continue
    edge_pairs.add(key)
    deduped.append(e)
print(f"  Unique STRING edges: {len(deduped)}")

edge_df = pd.DataFrame(deduped)
edge_df.to_csv(OUT / "04_string_edges.tsv", sep="\t", index=False)

# Build adjacency
adj = defaultdict(set)
for e in deduped:
    a, b = e["preferredName_A"], e["preferredName_B"]
    adj[a].add(b)
    adj[b].add(a)

# Per-set degree and shared-neighbor stats
sets = {k: data[k]["top"] for k in DATASETS}
deg_rows = []
nodes_in_graph = set(adj.keys())
for g in sorted(nodes_in_graph):
    sets_g = [k for k, s in sets.items() if g in s]
    nbrs = adj[g]
    nbr_by_set = {k: len(nbrs & s) for k, s in sets.items()}
    deg_rows.append({
        "gene": g,
        "in_sets": ",".join(sets_g),
        "degree": len(nbrs),
        **{f"nbrs_in_{k}": v for k, v in nbr_by_set.items()},
    })
deg_df = pd.DataFrame(deg_rows).sort_values("degree", ascending=False)
deg_df.to_csv(OUT / "05_string_node_degree.tsv", sep="\t", index=False)

# Cross-set "bridge" count: edges spanning two different sets
bridge_rows = []
for e in deduped:
    a, b = e["preferredName_A"], e["preferredName_B"]
    sets_a = [k for k, s in sets.items() if a in s]
    sets_b = [k for k, s in sets.items() if b in s]
    for sa in sets_a:
        for sb in sets_b:
            if sa < sb:
                bridge_rows.append({"set_A": sa, "set_B": sb, "gene_A": a, "gene_B": b, "score": e.get("score", "")})
            elif sa == sb:
                bridge_rows.append({"set_A": sa, "set_B": sa, "gene_A": a, "gene_B": b, "score": e.get("score", "")})
bridge_df = pd.DataFrame(bridge_rows)
bridge_df.to_csv(OUT / "06_string_bridge_edges.tsv", sep="\t", index=False)

# Connected components
def components(adj):
    seen = set()
    comps = []
    for start in adj:
        if start in seen:
            continue
        stack = [start]
        comp = []
        while stack:
            n = stack.pop()
            if n in seen:
                continue
            seen.add(n)
            comp.append(n)
            stack.extend(adj[n] - seen)
        comps.append(sorted(comp))
    return sorted(comps, key=len, reverse=True)

comps = components(adj)
print(f"  Connected components: {len(comps)}  (largest: {len(comps[0]) if comps else 0} nodes)")

comp_rows = []
for i, c in enumerate(comps):
    set_counts = Counter()
    for g in c:
        for k, s in sets.items():
            if g in s:
                set_counts[k] += 1
    comp_rows.append({
        "component_id": i, "size": len(c),
        **{f"n_{k}": set_counts.get(k, 0) for k in DATASETS},
        "genes": ",".join(c[:40]) + ("..." if len(c) > 40 else "")
    })
pd.DataFrame(comp_rows).to_csv(OUT / "07_string_components.tsv", sep="\t", index=False)

# Cross-set bridge summary
bridge_summary = bridge_df.groupby(["set_A", "set_B"]).size().reset_index(name="n_edges")
print("\n  Cross-set STRING edge counts:")
print(bridge_summary.to_string(index=False))
bridge_summary.to_csv(OUT / "08_string_bridge_summary.tsv", sep="\t", index=False)


# ============================================================
# 6. Top hub nodes per set
# ============================================================
print("\n  Top 10 hub nodes (by degree in STRING subgraph):")
print(deg_df.head(10).to_string(index=False))


# ============================================================
# Done
# ============================================================
print("\n" + "=" * 60)
print(f"All outputs under: {OUT}")
for f in sorted(OUT.glob("*.tsv")):
    n = sum(1 for _ in open(f)) - 1
    print(f"  {f.name}  ({n:,} rows)")
print("=" * 60)
