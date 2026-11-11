#!/bin/bash
# Submit from the repository root. The log paths below are relative, because
# Slurm does not expand variables in #SBATCH lines, and the folder must exist:
#   mkdir -p results/02_geneformer_in_silico_screen/logs
#   sbatch 02_geneformer_in_silico_screen/00_environment/00c_install_geneformer_pkg.sh
# Site-specific: partition, GPU constraint, node exclusion, array range, CPUs, memory and wall time are those of the original run (URI Unity cluster); adapt them to your scheduler.
#SBATCH --job-name=gf_pkg
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=results/02_geneformer_in_silico_screen/logs/gf_pkg_%j.out
#SBATCH --error=results/02_geneformer_in_silico_screen/logs/gf_pkg_%j.err
#SBATCH --mail-type=END,FAIL
# No --mail-user here: pass your own address on the command line (sbatch --mail-user=<address> ...).

# Install the Geneformer Python package (not model weights) by cloning the
# HuggingFace repo with GIT_LFS_SKIP_SMUDGE=1 (skips large binary/lfs files)
# then pip-installing only the Python source.
#
# PyPI deps (transformers, datasets, huggingface_hub, etc.) are already installed
# by the previous 00_install_and_download.sh job.
#
# Model weights are downloaded separately by 00d_download_model_weights.sh.

set -euo pipefail

# Repository and section configuration (gf_env.sh documents every GF_* variable).
# Under sbatch this script runs from a spool copy, so the repository is taken
# from the submit directory, or from NANOAI_REPO if that is exported.
GF_REPO="${NANOAI_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
if [ ! -f "${GF_REPO}/02_geneformer_in_silico_screen/gf_env.sh" ]; then
    echo "ERROR: repository not found at ${GF_REPO}; submit from the repository root or export NANOAI_REPO." >&2
    exit 1
fi
source "${GF_REPO}/02_geneformer_in_silico_screen/gf_env.sh"

PYTHON="${GF_PYTHON}"
PIP="${GF_PYTHON} -m pip"
CLONE_DIR="${GF_WORK}/models/geneformer_src"

echo "Job: ${SLURM_JOB_ID}  |  Node: ${SLURMD_NODENAME}  |  $(date)"

echo ""
echo "--- Cloning Geneformer (skip LFS binaries, Python source only) ---"
# GIT_LFS_SKIP_SMUDGE=1: clones text/source files only, skips .bin/.safetensors
rm -rf "${CLONE_DIR}"
GIT_LFS_SKIP_SMUDGE=1 git clone \
    https://huggingface.co/ctheodoris/Geneformer \
    "${CLONE_DIR}"

echo ""
echo "--- Cloned files ---"
ls "${CLONE_DIR}"

echo ""
echo "--- Installing Geneformer Python package ---"
${PIP} install --no-deps "${CLONE_DIR}"

echo ""
echo "--- Verifying geneformer import ---"
${PYTHON} -c "
import geneformer
from geneformer import TranscriptomeTokenizer, InSilicoPerturber, InSilicoPerturberStats
print('geneformer location:', geneformer.__file__)
print('TranscriptomeTokenizer OK')
print('InSilicoPerturber    OK')
print('InSilicoPerturberStats OK')
"

echo ""
echo "Done at $(date)"
echo "Next: sbatch 02_geneformer_in_silico_screen/00_environment/00e_fix_transformers_deps.sh, then 00d_download_model_weights.sh"
