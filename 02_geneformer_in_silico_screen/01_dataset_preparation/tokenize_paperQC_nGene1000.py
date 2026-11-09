#!/usr/bin/env python
"""Tokenize the paper-matched (n_genes>=1000) cell-type subsets for Geneformer V2."""
from pathlib import Path
from geneformer import TranscriptomeTokenizer

BASE = Path("/mnt/d/TBI_SC/Geneformer_TBI")
INPUTS = BASE / "h5ad/geneformer_inputs_paperQC_nGene1000"
OUTDIR = BASE / "tokenized_paperQC_nGene1000"
OUTDIR.mkdir(parents=True, exist_ok=True)

custom_attrs = {
    "condition": "condition",
    "cell_type": "cell_type",
    "sample": "sample",
    "donor": "donor",
    "region": "region",
    "age": "age",
    "sex": "sex",
    "n_counts": "n_counts",
}
PREFIXES = {
    "excitatory_neuron":  "GSE209552_excitatory_neurons_nGene1000",
    "inhibitory_neuron":  "GSE209552_inhibitory_neurons_nGene1000",
    "oligodendrocyte":    "GSE209552_oligodendrocytes_nGene1000",
    "OPC":                "GSE209552_OPCs_nGene1000",
}
for sub, prefix in PREFIXES.items():
    data_dir = INPUTS / sub
    out_dir = OUTDIR / sub
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== Tokenizing {sub} ===", flush=True)
    tk = TranscriptomeTokenizer(custom_attr_name_dict=custom_attrs, nproc=4, model_version="V2")
    tk.tokenize_data(
        data_directory=str(data_dir),
        output_directory=str(out_dir),
        output_prefix=prefix,
        file_format="h5ad",
    )
    print(f"  wrote {out_dir}/{prefix}.dataset", flush=True)
print("\nAll paper-matched subsets tokenized.")
