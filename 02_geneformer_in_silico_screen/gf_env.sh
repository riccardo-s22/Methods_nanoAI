#!/usr/bin/env bash
# Path configuration for Results section 02 (Geneformer in-silico screen), for
# shell and Slurm scripts. Source it, do not execute it:
#
#     source 02_geneformer_in_silico_screen/gf_env.sh
#
# It is the shell twin of gf_paths.py: both resolve the same variables to the
# same defaults, so a value exported here (or before sourcing) is what the
# Python scripts see. Built on shared/env.sh (NANOAI_REPO, NANOAI_DATA,
# NANOAI_RESULTS). Nothing here is an analysis parameter: only locations.
#
# Every variable is respected if already set; otherwise it gets its default.
#
#   GF_WORK            Project working root of the FT (single-cell-type) pipeline.
#                      Default: <repo>/work/Geneformer_TBI   (work/ is git-ignored)
#                      The layout below it is the layout of the original project
#                      tree, unchanged, so pointing GF_WORK at a copy of that
#                      tree reproduces the original behaviour:
#                          h5ad/                          merged, QC'd and per-cell-type AnnData
#                          metadata/  reports/            sample metadata, QC reports, sample manifest
#                          tokenized_paperQC_nGene1000/   rank-value encoded .dataset folders
#                          models/Geneformer/             pretrained weights + fine-tuned checkpoints
#                          results/                       fine-tuning, state embeddings, ISP chunks, statistics
#   GF_JOINTPG         Project root of the MT (multi-task, "JointPG") model, whose
#                      result tables the MT comparisons read.
#                      Default: <repo>/work/Geneformer_JointPG
#   GF_MECHANO         Root of the MechanoGEPred project; only
#                      data/mechano_genes_list.txt is read from it.
#                      Default: <repo>/work/MechanoGEPred
#   GF_GEO             Where the GEO GSE209552 download lives.
#                      Default: <repo>/work/GSE209552
#                          GSE209552_RAW/          extracted GSE209552_RAW.tar
#                          GSE209552_family.soft   series SOFT file (or a folder of
#                                                  that name holding the file)
#                          SraRunTable.csv         SRA Run Selector metadata table
#   GF_SOFT            The SOFT file itself.          Default: derived from GF_GEO
#   GF_SRA_RUN_TABLE   The SRA run table itself.      Default: derived from GF_GEO
#   GF_MODEL_DIR       Snapshot of the Hugging Face repository ctheodoris/Geneformer
#                      (weights, token dictionaries, package source).
#                      Default: $GF_WORK/models/Geneformer
#   GF_PRETRAINED      Pretrained backbone that is fine-tuned.
#                      Default: $GF_MODEL_DIR/Geneformer-V2-104M
#   GF_GENEFORMER_PKG  Folder holding token_dictionary_gc104M.pkl and
#                      gene_name_id_dict_gc104M.pkl.
#                      Default: $GF_MODEL_DIR/geneformer
#   GF_FT_MODELS       Parent folder Geneformer writes the fine-tuned FT classifier into.
#                      Default: $GF_MODEL_DIR/fine_tuned_models_tbi
#   GF_FT_MODEL        A fine-tuned FT checkpoint to use directly. When set it
#                      replaces the path stored in
#                      $GF_WORK/results/finetune_tbi_excitatory/trained_model_path.txt.
#                      Default: unset (the pointer file is used)
#   GF_MT_MODEL        The retained MT checkpoint (published: step 600). Not available.
#                      Default: $GF_JOINTPG/models/garza_mt/ksplit1/checkpoint-600
#   GF_CONDA_ENV       Prefix of the conda environment holding Geneformer and PyTorch.
#                      Default: unset
#   GF_PYTHON          Interpreter for everything that imports geneformer or torch.
#                      Default: $GF_CONDA_ENV/bin/python if GF_CONDA_ENV is set,
#                      otherwise the `python` on PATH (the activated environment).
#   GF_PREP_PYTHON     Interpreter for the scanpy/harmonypy/scrublet preparation
#                      scripts, which were run in a separate environment.
#                      Default: $GF_PYTHON
#   GF_LOGS            Log folder for non-Slurm runs.
#                      Default: $NANOAI_RESULTS/02_geneformer_in_silico_screen/logs
#
# Also exported, derived and not meant to be overridden:
#   GF_SECTION_DIR     <repo>/02_geneformer_in_silico_screen
#   GF_DATA            $NANOAI_DATA/02_geneformer_in_silico_screen
#   GF_RESULTS         $NANOAI_RESULTS/02_geneformer_in_silico_screen
#
# Slurm: `#SBATCH --output/--error` cannot expand variables, so the Slurm scripts
# of this section use the relative path
# results/02_geneformer_in_silico_screen/logs/ and are submitted from the
# repository root. The folder has to exist before submission:
#
#     mkdir -p results/02_geneformer_in_silico_screen/logs
#     sbatch 02_geneformer_in_silico_screen/<stage>/<script>.sh
#
# Under sbatch a script runs from a spool copy and cannot locate itself, so the
# Slurm scripts find the repository through SLURM_SUBMIT_DIR (or NANOAI_REPO, if
# exported) and then source this file.

GF_SECTION_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export GF_SECTION_DIR
# shellcheck source=../shared/env.sh
source "${GF_SECTION_DIR}/../shared/env.sh"

_gf_abs() {  # absolute path without requiring the target to exist
    case "$1" in
        /*) printf '%s\n' "$1" ;;
        *)  printf '%s\n' "$(pwd)/$1" ;;
    esac
}

GF_WORK="$(_gf_abs "${GF_WORK:-${NANOAI_REPO}/work/Geneformer_TBI}")"
GF_JOINTPG="$(_gf_abs "${GF_JOINTPG:-${NANOAI_REPO}/work/Geneformer_JointPG}")"
GF_MECHANO="$(_gf_abs "${GF_MECHANO:-${NANOAI_REPO}/work/MechanoGEPred}")"
GF_GEO="$(_gf_abs "${GF_GEO:-${NANOAI_REPO}/work/GSE209552}")"
export GF_WORK GF_JOINTPG GF_MECHANO GF_GEO

GF_MODEL_DIR="$(_gf_abs "${GF_MODEL_DIR:-${GF_WORK}/models/Geneformer}")"
GF_PRETRAINED="$(_gf_abs "${GF_PRETRAINED:-${GF_MODEL_DIR}/Geneformer-V2-104M}")"
GF_GENEFORMER_PKG="$(_gf_abs "${GF_GENEFORMER_PKG:-${GF_MODEL_DIR}/geneformer}")"
GF_FT_MODELS="$(_gf_abs "${GF_FT_MODELS:-${GF_MODEL_DIR}/fine_tuned_models_tbi}")"
GF_MT_MODEL="$(_gf_abs "${GF_MT_MODEL:-${GF_JOINTPG}/models/garza_mt/ksplit1/checkpoint-600}")"
export GF_MODEL_DIR GF_PRETRAINED GF_GENEFORMER_PKG GF_FT_MODELS GF_MT_MODEL
if [ -n "${GF_FT_MODEL:-}" ]; then
    GF_FT_MODEL="$(_gf_abs "${GF_FT_MODEL}")"
    export GF_FT_MODEL
fi

# GEO inputs. The workstation copy held the SOFT file inside a folder of the same
# name and the run table under a browser-duplicate name; accept both forms.
if [ -z "${GF_SOFT:-}" ]; then
    if [ -d "${GF_GEO}/GSE209552_family.soft" ]; then
        GF_SOFT="${GF_GEO}/GSE209552_family.soft/GSE209552_family.soft"
    else
        GF_SOFT="${GF_GEO}/GSE209552_family.soft"
    fi
fi
if [ -z "${GF_SRA_RUN_TABLE:-}" ]; then
    if [ ! -e "${GF_GEO}/SraRunTable.csv" ] && [ -e "${GF_GEO}/SraRunTable (1).csv" ]; then
        GF_SRA_RUN_TABLE="${GF_GEO}/SraRunTable (1).csv"
    else
        GF_SRA_RUN_TABLE="${GF_GEO}/SraRunTable.csv"
    fi
fi
GF_SOFT="$(_gf_abs "${GF_SOFT}")"
GF_SRA_RUN_TABLE="$(_gf_abs "${GF_SRA_RUN_TABLE}")"
export GF_SOFT GF_SRA_RUN_TABLE

# Interpreters.
export GF_CONDA_ENV="${GF_CONDA_ENV:-}"
if [ -z "${GF_PYTHON:-}" ]; then
    if [ -n "${GF_CONDA_ENV}" ]; then
        GF_PYTHON="${GF_CONDA_ENV}/bin/python"
    else
        GF_PYTHON="python"
    fi
fi
export GF_PYTHON
export GF_PREP_PYTHON="${GF_PREP_PYTHON:-${GF_PYTHON}}"

# Shipped data, rerun outputs, logs.
export GF_DATA="${NANOAI_DATA}/02_geneformer_in_silico_screen"
export GF_RESULTS="${NANOAI_RESULTS}/02_geneformer_in_silico_screen"
GF_LOGS="$(_gf_abs "${GF_LOGS:-${GF_RESULTS}/logs}")"
export GF_LOGS

unset -f _gf_abs
