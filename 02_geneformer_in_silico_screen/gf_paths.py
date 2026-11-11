"""Path configuration for Results section 02 (Geneformer in-silico screen).

This file and gf_env.sh (its shell twin) are the single place where the section
is pointed at a machine. Nothing here is an analysis parameter: only locations.
Every variable is read from the environment and falls back to a default inside
the repository, so a fresh clone needs no editing.

Built on shared/repo_paths.py (NANOAI_DATA -> <repo>/data, NANOAI_RESULTS ->
<repo>/results), which it re-exports.

-----------------------------------------------------------------------------
Variables (all optional; set them in the environment before running)
-----------------------------------------------------------------------------
GF_WORK            Project working root of the FT (single-cell-type) pipeline.
                   Default: <repo>/work/Geneformer_TBI   (work/ is git-ignored)
                   The layout below it is the layout of the original project
                   tree, unchanged, so pointing GF_WORK at a copy of that tree
                   reproduces the original behaviour:
                       h5ad/                          merged, QC'd and per-cell-type AnnData
                       metadata/  reports/            sample metadata, QC reports, sample manifest
                       tokenized_paperQC_nGene1000/   rank-value encoded .dataset folders
                       models/Geneformer/             pretrained weights + fine-tuned checkpoints
                       results/                       fine-tuning, state embeddings, ISP chunks, statistics
GF_JOINTPG         Project root of the MT (multi-task, "JointPG") model, whose
                   result tables the MT comparisons read.
                   Default: <repo>/work/Geneformer_JointPG
GF_MECHANO         Root of the MechanoGEPred project; only
                   data/mechano_genes_list.txt is read from it.
                   Default: <repo>/work/MechanoGEPred
GF_GEO             Where the GEO GSE209552 download lives.
                   Default: <repo>/work/GSE209552
                       GSE209552_RAW/            extracted GSE209552_RAW.tar
                       GSE209552_family.soft     series SOFT file (or a folder of
                                                 that name holding the file)
                       SraRunTable.csv           SRA Run Selector metadata table
GF_SOFT            The SOFT file itself.            Default: derived from GF_GEO
GF_SRA_RUN_TABLE   The SRA run table itself.        Default: derived from GF_GEO
GF_MODEL_DIR       Snapshot of the Hugging Face repository ctheodoris/Geneformer
                   (weights, token dictionaries, package source).
                   Default: $GF_WORK/models/Geneformer
GF_PRETRAINED      Pretrained backbone that is fine-tuned.
                   Default: $GF_MODEL_DIR/Geneformer-V2-104M
GF_GENEFORMER_PKG  Folder holding token_dictionary_gc104M.pkl and
                   gene_name_id_dict_gc104M.pkl.
                   Default: $GF_MODEL_DIR/geneformer
GF_FT_MODELS       Parent folder Geneformer writes the fine-tuned FT classifier into.
                   Default: $GF_MODEL_DIR/fine_tuned_models_tbi
GF_FT_MODEL        A fine-tuned FT checkpoint to use directly. When set it
                   replaces the path stored in
                   $GF_WORK/results/finetune_tbi_excitatory/trained_model_path.txt.
                   Default: unset (the pointer file is used)
GF_MT_MODEL        The retained MT checkpoint (published: step 600). Not
                   available; see MT_RECON below.
                   Default: $GF_JOINTPG/models/garza_mt/ksplit1/checkpoint-600
GF_CONDA_ENV       Prefix of the conda environment holding Geneformer and PyTorch.
                   Default: unset
GF_PYTHON          Interpreter for everything that imports geneformer or torch.
                   Default: $GF_CONDA_ENV/bin/python if GF_CONDA_ENV is set,
                   otherwise the `python` on PATH (the activated environment).
GF_PREP_PYTHON     Interpreter for the scanpy/harmonypy/scrublet preparation
                   scripts, which were run in a separate environment.
                   Default: $GF_PYTHON
GF_LOGS            Log folder for non-Slurm runs.
                   Default: $NANOAI_RESULTS/02_geneformer_in_silico_screen/logs
                   (Slurm scripts write to the same folder through a relative
                   path, which is why they are submitted from the repository root.)

GF_PYTHON, GF_PREP_PYTHON and GF_CONDA_ENV matter to the shell scripts only;
they are mirrored here so both files document the same set.

Where things are written. The GPU-side pipeline (stages 00 to 03, and stages 07
and 09, which read and write the same per-arm folders) works inside GF_WORK, as
it did on the cluster. The CPU analyses after it (04_pathway_enrichment except
09, 05_drug_repurposing, 06_confounders_and_robustness, figures) write below
<results>/02_geneformer_in_silico_screen/ and find their inputs through
analysis_input() at the end of this file.

-----------------------------------------------------------------------------
Usage, from a script N directories below this section directory
-----------------------------------------------------------------------------
    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    *[".."] * N))
    import gf_paths
    BASE = Path(gf_paths.GF_WORK)

Importing this module creates no directory and touches no file.
"""
import os
import sys

SECTION_DIR = os.path.dirname(os.path.abspath(__file__))
SECTION = os.path.basename(SECTION_DIR)  # "02_geneformer_in_silico_screen"

sys.path.insert(0, os.path.join(SECTION_DIR, "..", "shared"))
from repo_paths import (  # noqa: E402,F401  (re-exported for the section's scripts)
    DATA_DIR,
    REPO_ROOT,
    RESULTS_DIR,
    data_path,
    env_path,
    repo_path,
    results_dir,
)


def _root(name, default):
    """Absolute path from environment variable `name`, else `default`."""
    return os.path.abspath(os.environ.get(name) or default)


# ---------------------------------------------------------------------------
# Roots
# ---------------------------------------------------------------------------
GF_WORK = _root("GF_WORK", repo_path("work", "Geneformer_TBI"))
GF_JOINTPG = _root("GF_JOINTPG", repo_path("work", "Geneformer_JointPG"))
GF_MECHANO = _root("GF_MECHANO", repo_path("work", "MechanoGEPred"))
GF_GEO = _root("GF_GEO", repo_path("work", "GSE209552"))

GF_MODEL_DIR = _root("GF_MODEL_DIR", os.path.join(GF_WORK, "models", "Geneformer"))
GF_PRETRAINED = _root("GF_PRETRAINED", os.path.join(GF_MODEL_DIR, "Geneformer-V2-104M"))
GF_GENEFORMER_PKG = _root("GF_GENEFORMER_PKG", os.path.join(GF_MODEL_DIR, "geneformer"))
GF_FT_MODELS = _root("GF_FT_MODELS", os.path.join(GF_MODEL_DIR, "fine_tuned_models_tbi"))
GF_FT_MODEL = os.path.abspath(os.environ["GF_FT_MODEL"]) if os.environ.get("GF_FT_MODEL") else ""
GF_MT_MODEL = _root(
    "GF_MT_MODEL", os.path.join(GF_JOINTPG, "models", "garza_mt", "ksplit1", "checkpoint-600")
)

GF_CONDA_ENV = os.environ.get("GF_CONDA_ENV") or ""
GF_PYTHON = os.environ.get("GF_PYTHON") or (
    os.path.join(GF_CONDA_ENV, "bin", "python") if GF_CONDA_ENV else "python"
)
GF_PREP_PYTHON = os.environ.get("GF_PREP_PYTHON") or GF_PYTHON

# Shipped source data and the outputs of the CPU-side reruns for this section.
GF_DATA = data_path(SECTION)
GF_RESULTS = os.path.join(RESULTS_DIR, SECTION)  # not created here; use gf_results_dir()
GF_LOGS = _root("GF_LOGS", os.path.join(GF_RESULTS, "logs"))

# Where the section's importable helper modules live.
FT_FINETUNE_CODE = os.path.join(SECTION_DIR, "02_model_finetuning", "ft_excitatory")
FT_ISP_CODE = os.path.join(SECTION_DIR, "03_in_silico_perturbation", "ft_excitatory")


def work_path(*parts):
    """Path below GF_WORK. Not created, not checked."""
    return os.path.join(GF_WORK, *parts)


def jointpg_path(*parts):
    """Path below GF_JOINTPG. Not created, not checked."""
    return os.path.join(GF_JOINTPG, *parts)


def gf_results_dir(*parts):
    """Output directory below <results>/02_geneformer_in_silico_screen, created."""
    return results_dir(SECTION, *parts)


def first_existing(*candidates):
    """The first candidate that exists; the first candidate if none does, so an
    error message names the primary location."""
    for path in candidates:
        if os.path.exists(path):
            return path
    return candidates[0]


def add_ft_code_to_syspath():
    """Make finetuned_celltype_config, perturb_safeguards and
    biological_gene_classes importable. The scripts were written for one flat
    scripts/ folder; in this repository those modules sit in two folders."""
    for folder in (FT_ISP_CODE, FT_FINETUNE_CODE):
        if folder not in sys.path:
            sys.path.insert(0, folder)


# ---------------------------------------------------------------------------
# GEO download (inputs of 01_dataset_preparation)
# ---------------------------------------------------------------------------
GEO_RAW = os.path.join(GF_GEO, "GSE209552_RAW")


def _geo_soft():
    # The workstation copy was a folder GSE209552_family.soft/ holding the file;
    # `gunzip GSE209552_family.soft.gz` gives the file directly. Accept both.
    soft = os.path.join(GF_GEO, "GSE209552_family.soft")
    return os.path.join(soft, "GSE209552_family.soft") if os.path.isdir(soft) else soft


def _geo_sra_table():
    # The workstation copy carried a browser-duplicate name.
    return first_existing(
        os.path.join(GF_GEO, "SraRunTable.csv"),
        os.path.join(GF_GEO, "SraRunTable (1).csv"),
    )


GF_SOFT = _root("GF_SOFT", _geo_soft())
GF_SRA_RUN_TABLE = _root("GF_SRA_RUN_TABLE", _geo_sra_table())


# ---------------------------------------------------------------------------
# Token dictionaries
# ---------------------------------------------------------------------------
TOKEN_DICT = os.path.join(GF_GENEFORMER_PKG, "token_dictionary_gc104M.pkl")
GENE_NAME_DICT = os.path.join(GF_GENEFORMER_PKG, "gene_name_id_dict_gc104M.pkl")


# ---------------------------------------------------------------------------
# FT model: fine-tuned checkpoint and in-silico perturbation tables
# ---------------------------------------------------------------------------
FT_PREP_DIR = os.path.join(GF_WORK, "results", "finetune_tbi_excitatory")
FT_POINTER = os.path.join(FT_PREP_DIR, "trained_model_path.txt")
# The checkpoint of the published run (datestamped by Geneformer on 2026-05-15).
FT_PUBLISHED_CHECKPOINT = os.path.join(
    GF_FT_MODELS, "260515_geneformer_cellClassifier_excit_tbi_clf", "ksplit1", "checkpoint-306"
)


def resolve_model_pointer(pointer_file=FT_POINTER):
    """Directory of the fine-tuned FT checkpoint.

    GF_FT_MODEL wins when set. Otherwise the path stored in the pointer file is
    used. That file records an absolute path of the machine the model was
    trained on; if that path does not exist here, the part below
    models/Geneformer/ is looked up under GF_MODEL_DIR, and the part below the
    project root under GF_WORK, so a copied project tree keeps working. If
    nothing exists the stored path is returned unchanged and the caller reports
    it missing.
    """
    if GF_FT_MODEL:
        return GF_FT_MODEL
    with open(pointer_file) as handle:
        stored = handle.read().strip()
    if os.path.exists(stored):
        return stored
    norm = stored.replace("\\", "/")
    for marker, root in (("/models/Geneformer/", GF_MODEL_DIR), ("/Geneformer_TBI/", GF_WORK)):
        if marker in norm:
            candidate = os.path.join(root, *norm.split(marker, 1)[1].split("/"))
            if os.path.exists(candidate):
                return candidate
    return stored


# arm -> (run folder below $GF_WORK/results, perturbation mode, table prefix)
FT_ARMS = {
    "rescue_delete": (
        "perturb_neurons_finetuned_tbi_to_control", "delete", "excit_finetuned_rescue_delete"),
    "rescue_overexpress": (
        "perturb_neurons_finetuned_tbi_to_control", "overexpress", "excit_finetuned_rescue_overexpress"),
    "induction_delete": (
        "perturb_neurons_finetuned_control_to_tbi", "delete", "excit_finetuned_delete"),
    "induction_overexpress": (
        "perturb_neurons_finetuned_control_to_tbi", "overexpress", "excit_finetuned_overexpress"),
}
FT_STATE_EMBS = os.path.join(
    GF_WORK, "results", "perturb_neurons_finetuned_control_to_tbi", "state_embs",
    "excitatory_finetuned_state_embs.pkl",
)


def ft_run_root(arm):
    return os.path.join(GF_WORK, "results", FT_ARMS[arm][0])


def ft_stats_work(arm, kind="unfiltered_all_gene_ranked"):
    """Per-gene table of one FT arm as the pipeline writes it below GF_WORK."""
    run, mode, prefix = FT_ARMS[arm]
    return os.path.join(GF_WORK, "results", run, f"stats_{mode}", f"{prefix}_{kind}.tsv")


def ft_stats_shipped(arm, kind="unfiltered_all_gene_ranked"):
    """The same table as shipped in data/ (unfiltered_all_gene_ranked and
    filtered_high_confidence_gene_ranked only)."""
    return os.path.join(
        GF_DATA, "in_silico_perturbation", "ft_excitatory", arm, f"{FT_ARMS[arm][2]}_{kind}.tsv"
    )


def ft_stats_table(arm, kind="unfiltered_all_gene_ranked"):
    """A fresh rerun below GF_WORK if there is one, otherwise the shipped table."""
    return first_existing(ft_stats_work(arm, kind), ft_stats_shipped(arm, kind))


FT_SHIPPED_KINDS = ("unfiltered_all_gene_ranked", "filtered_high_confidence_gene_ranked")


def seed_ft_stats_from_shipped():
    """Copy the shipped per-gene tables of the four FT arms from data/ into the
    GF_WORK layout (results/<run>/stats_<mode>/), where stage 07 and the scripts
    that read the project tree expect them. This lets everything after the GPU
    stages run from a fresh clone. data/ is only read; a table already present
    below GF_WORK (a real rerun) is never overwritten. Returns the copied paths.
    """
    import shutil

    copied = []
    for arm in FT_ARMS:
        for kind in FT_SHIPPED_KINDS:
            src, dst = ft_stats_shipped(arm, kind), ft_stats_work(arm, kind)
            if os.path.exists(dst) or not os.path.exists(src):
                continue
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copyfile(src, dst)
            copied.append(dst)
    return copied


# ---------------------------------------------------------------------------
# MT model: tables of the published run that downstream scripts read
# ---------------------------------------------------------------------------
# These follow the layout of the original MT project. None of them is shipped:
# the MT in-silico perturbation tables are not available.
MT_ISP = os.path.join(GF_JOINTPG, "results", "isp")
MT_DIRECTIONS = {"rescue": "TBI_to_control", "induction": "control_to_TBI"}
MT_DETECTION_TABLE = os.path.join(
    MT_ISP, "garza_disease_binary_TBI_to_control", "mechano_analysis", "layer3",
    "garza_qc4_de_by_celltype.tsv",
)
MECHANO_GENES = os.path.join(GF_MECHANO, "data", "mechano_genes_list.txt")


def mt_stats_table(direction, mode, cell_type="excitatory_neuron",
                   kind="unfiltered_all_gene_ranked"):
    """Per-cell-type table of the published MT screen.
    direction: "TBI_to_control" or "control_to_TBI"; mode: "delete" or "overexpress"."""
    run = f"garza_disease_binary_{direction}"
    return os.path.join(
        MT_ISP, run, f"stats_{mode}", "by_cell_type", f"isp_{run}_{mode}_{cell_type}_{kind}.tsv"
    )


# ---------------------------------------------------------------------------
# MT model: default locations for the *_RECONSTRUCTED scripts
# ---------------------------------------------------------------------------
# The reconstructed MT scripts take every location on the command line. The
# defaults below let them run without arguments. They are a layout chosen for
# this repository, NOT the layout of the original MT project, which was not
# recorded apart from the tables listed above.
MT_RECON = {
    # mt_01 input: the paper-matched primary object written by 01_dataset_preparation
    "primary_h5ad": os.path.join(
        GF_WORK, "h5ad", "Garza_GSE209552_all_raw_counts_paperQC_nGene1000_annotated.h5ad"),
    "depth_matched_h5ad": os.path.join(GF_JOINTPG, "h5ad", "garza_mt_depth_matched.h5ad"),
    "split_h5ad": os.path.join(GF_JOINTPG, "h5ad", "garza_mt_depth_matched_split.h5ad"),
    # Tokenized corpora. No script in the repository writes these.
    "tokenized_train": os.path.join(GF_JOINTPG, "tokenized", "garza_mt_train.dataset"),
    "tokenized_valid": os.path.join(GF_JOINTPG, "tokenized", "garza_mt_valid.dataset"),
    "tokenized_test": os.path.join(GF_JOINTPG, "tokenized", "garza_mt_test.dataset"),
    "tokenized_all": os.path.join(GF_JOINTPG, "tokenized", "garza_mt_all.dataset"),
    "model_save": os.path.join(GF_JOINTPG, "models", "garza_mt"),
    "finetune_results": os.path.join(GF_JOINTPG, "results", "finetune_garza_mt"),
    "state_embs_dir": os.path.join(GF_JOINTPG, "results", "isp_RECONSTRUCTED", "state_embs"),
    "state_embs": os.path.join(
        GF_JOINTPG, "results", "isp_RECONSTRUCTED", "state_embs", "garza_mt_state_embs.pkl"),
    "isp_chunks": os.path.join(GF_JOINTPG, "results", "isp_RECONSTRUCTED", "chunks"),
    "isp_stats": os.path.join(GF_JOINTPG, "results", "isp_RECONSTRUCTED", "stats"),
    "umap": os.path.join(GF_JOINTPG, "results", "isp_RECONSTRUCTED", "umap"),
}


# ---------------------------------------------------------------------------
# Downstream CPU analyses: 04_pathway_enrichment, 05_drug_repurposing,
# 06_confounders_and_robustness, figures
# ---------------------------------------------------------------------------
# These scripts write below <results>/02_geneformer_in_silico_screen/, in the
# folder names the original project used below its results/ folder
# (pathway_heatmap/, drug_repurpose/, cross_direction_network_v2/, umap/), plus
# confounders/ and figures/ for files that had no folder of their own.
#
# What one of them reads from an earlier one is looked up in three places, in
# this order, so that a rerun feeds the next step, a copy of the original
# project tree behaves as it did, and a fresh clone falls back on data/:
#   1. <results>/02_geneformer_in_silico_screen/<parts>   written by a rerun here
#   2. $GF_WORK/results/<parts>                            the original project tree
#   3. the table shipped in data/, where one exists
SHIPPED_PATHWAY_MATRIX = {
    model: os.path.join(GF_DATA, "pathway_enrichment", f"excit_pathway_enrichment_matrix_{model}.tsv")
    for model in ("FT", "MT")
}
SHIPPED_DRUG_DIR = os.path.join(GF_DATA, "drug_signature_enrichment")


def analysis_input(*parts, shipped=None):
    """A file or folder written by an earlier CPU analysis of this section.
    Returns the first existing of the three locations listed above (the first
    one if none exists, so an error names where a rerun would put it)."""
    candidates = [os.path.join(GF_RESULTS, *parts), os.path.join(GF_WORK, "results", *parts)]
    if shipped:
        candidates.append(shipped)
    return first_existing(*candidates)


def drug_table(*parts):
    """A table of the drug-signature screen, e.g. drug_table("FT",
    "drug_candidates_ranked_union.tsv") or drug_table("crossmodel_consensus.tsv").
    The shipped copies keep the original file names below
    data/02_geneformer_in_silico_screen/drug_signature_enrichment/."""
    return analysis_input("drug_repurpose", *parts, shipped=os.path.join(SHIPPED_DRUG_DIR, *parts))


def shipped_pathway_arm(label):
    """Enrichr adjusted P-values of one arm (label such as "FT_delete" or
    "MT_oe_c2t") read from the shipped Supplementary Figure 6 source matrix, in
    the three columns the heatmap scripts read from <label>_enrichr_allterms.csv
    (Gene_set, Term, Adjusted P-value). Returns None when the matrix is not
    there or holds no value for that arm. Rows without a value for the arm are
    terms Enrichr did not return for it and are left out, as they are absent
    from the per-arm file. The FT matrix holds every term; the MT matrix holds
    the terms significant in at least one arm, which is all its heatmap draws.
    """
    import pandas as pd

    path = SHIPPED_PATHWAY_MATRIX.get(label.split("_", 1)[0])
    if not path or not os.path.exists(path):
        return None
    # round_trip: the values must be the exact doubles the matrix was written
    # from; the default parser can be off by one unit in the last place.
    matrix = pd.read_csv(path, sep="\t", float_precision="round_trip")
    if label not in matrix.columns:
        return None
    arm = matrix[["Gene_set", "Term", label]].dropna(subset=[label])
    if arm.empty:
        return None
    return arm.rename(columns={label: "Adjusted P-value"}).reset_index(drop=True)


if __name__ == "__main__":
    # `python gf_paths.py` prints the resolved configuration.
    # `python gf_paths.py --seed-ft-stats` first copies the shipped FT tables into GF_WORK.
    if "--seed-ft-stats" in sys.argv[1:]:
        for _path in seed_ft_stats_from_shipped():
            print("seeded", _path)
    for _name in (
        "REPO_ROOT", "DATA_DIR", "RESULTS_DIR", "GF_WORK", "GF_JOINTPG", "GF_MECHANO",
        "GF_GEO", "GF_SOFT", "GF_SRA_RUN_TABLE", "GF_MODEL_DIR", "GF_PRETRAINED",
        "GF_GENEFORMER_PKG", "GF_FT_MODELS", "GF_FT_MODEL", "GF_MT_MODEL", "GF_CONDA_ENV",
        "GF_PYTHON", "GF_PREP_PYTHON", "GF_DATA", "GF_RESULTS", "GF_LOGS",
    ):
        _value = globals()[_name]
        _mark = "" if not str(_value).startswith(os.sep) else (
            "" if os.path.exists(_value) else "   [missing]")
        print(f"{_name:18s} {_value or '(unset)'}{_mark}")
