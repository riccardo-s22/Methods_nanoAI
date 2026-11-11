#!/usr/bin/env python
"""
Count genes detected in excitatory and inhibitory neurons.
Reports:
  - number of cells per cell type / condition
  - number of genes detected (non-zero in >= 1 cell)
  - number of those genes present in the Geneformer token dictionary
  - planned chunk sizes for the full ISP array job

Run on a login node (lightweight, CPU only):
  python scripts/01_count_neuronal_genes.py
"""
import json
from pathlib import Path

import anndata as ad
import numpy as np

# Section path configuration (GF_WORK, GF_GEO, ...): see ../gf_paths.py.
import os, sys  # noqa: E401
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gf_paths  # noqa: E402

BASE       = Path(gf_paths.GF_WORK)
SUBSETS    = BASE / "h5ad/celltype_subsets"
TOKENIZED  = BASE / "tokenized"
MODEL_DIR  = Path(gf_paths.GF_MODEL_DIR)
N_CHUNKS   = 6   # desired number of Slurm array tasks for excitatory neurons

def count_detected_genes(h5ad_path: Path) -> dict:
    adata = ad.read_h5ad(h5ad_path)
    import scipy.sparse as sp
    X = adata.X
    if sp.issparse(X):
        gene_counts = np.array(X.sum(axis=0)).flatten()
    else:
        gene_counts = X.sum(axis=0)
    detected_mask = gene_counts > 0
    detected_ids  = adata.var.loc[detected_mask, "ensembl_id"].tolist()
    return {
        "n_cells":        adata.n_obs,
        "n_genes_total":  adata.n_vars,
        "n_genes_detected": int(detected_mask.sum()),
        "detected_ensembl_ids": detected_ids,
        "condition_counts": adata.obs["condition"].value_counts().to_dict(),
    }

def load_gf_token_dict(model_dir: Path) -> dict | None:
    """Load Geneformer gene -> token mapping for the V2-104M model."""
    import pickle
    # Prefer the gc104M token dictionary (matches our V2-104M model)
    preference = [
        model_dir / "geneformer" / "token_dictionary_gc104M.pkl",
        model_dir / "geneformer" / "gene_name_id_dict_gc104M.pkl",
    ]
    candidates = [p for p in preference if p.exists()]
    # fallback: any token dict
    if not candidates:
        candidates = (
            list(model_dir.rglob("token_dictionary*.pkl")) +
            list(model_dir.rglob("gene_name_id_dict*.pkl"))
        )
    if not candidates:
        print("  [WARN] Geneformer token dictionary not found.")
        print(f"         Expected in: {model_dir}/geneformer/")
        print("         Run 00d_download_model_weights.sh first.")
        return None
    with open(candidates[0], "rb") as f:
        d = pickle.load(f)
    print(f"  Token dict loaded from: {candidates[0].name}  ({len(d):,} entries)")
    return d

print("=" * 70)
print("Geneformer TBI — Neuronal Gene Detection Report")
print("=" * 70)

results = {}
for cell_type, fname in [
    ("excitatory_neuron",  "excitatory_neurons_raw_geneformer_ready.h5ad"),
    ("inhibitory_neuron",  "inhibitory_neurons_raw_geneformer_ready.h5ad"),
]:
    path = SUBSETS / fname
    if not path.exists():
        print(f"\n[SKIP] {fname} not found")
        continue
    print(f"\n--- {cell_type} ---")
    info = count_detected_genes(path)
    results[cell_type] = info
    print(f"  Cells:               {info['n_cells']:,}")
    print(f"  Condition counts:    {info['condition_counts']}")
    print(f"  Genes in matrix:     {info['n_genes_total']:,}")
    print(f"  Genes detected:      {info['n_genes_detected']:,}")

# Try to intersect with Geneformer vocabulary
print(f"\n--- Geneformer token dictionary ---")
token_dict = load_gf_token_dict(MODEL_DIR)

for cell_type, info in results.items():
    detected = set(info["detected_ensembl_ids"])
    if token_dict is not None:
        in_vocab = detected & set(token_dict.keys())
        n_vocab  = len(in_vocab)
        print(f"  {cell_type}: {n_vocab:,} / {len(detected):,} detected genes in GF vocabulary")
        info["n_genes_in_vocab"] = n_vocab
        info["genes_in_vocab"]   = sorted(in_vocab)
    else:
        info["n_genes_in_vocab"] = info["n_genes_detected"]
        info["genes_in_vocab"]   = info["detected_ensembl_ids"]

# Chunk plan for excitatory neurons
print(f"\n--- Chunk plan (excitatory neurons, N_CHUNKS={N_CHUNKS}) ---")
if "excitatory_neuron" in results:
    ex = results["excitatory_neuron"]
    n_genes = ex["n_genes_in_vocab"]
    chunk_size = int(np.ceil(n_genes / N_CHUNKS))
    print(f"  Genes to perturb:  {n_genes:,}")
    print(f"  N chunks:          {N_CHUNKS}")
    print(f"  Genes per chunk:   ~{chunk_size:,}")
    print(f"  Slurm array:       --array=0-{N_CHUNKS - 1}")

    # Save gene list for chunks
    gene_list_path = BASE / "results/perturb_neurons/excitatory_neuron_genes_to_perturb.json"
    gene_list_path.parent.mkdir(parents=True, exist_ok=True)  # fresh GF_WORK: folder may not exist yet
    with open(gene_list_path, "w") as f:
        json.dump({
            "cell_type": "excitatory_neuron",
            "n_genes": n_genes,
            "n_chunks": N_CHUNKS,
            "chunk_size": chunk_size,
            "genes": ex["genes_in_vocab"],
        }, f, indent=2)
    print(f"\n  Gene list saved to: {gene_list_path}")

print(f"\n  Cells (excitatory):  {results.get('excitatory_neuron', {}).get('n_cells', 'N/A'):,}")
print(f"  Cells (inhibitory):  {results.get('inhibitory_neuron', {}).get('n_cells', 'N/A'):,}")
print()
