#!/bin/bash
#SBATCH --job-name=gf_pkg
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/gf_pkg_%j.out
#SBATCH --error=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/gf_pkg_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=riccardo-s@uri.edu

# Install the Geneformer Python package (not model weights) by cloning the
# HuggingFace repo with GIT_LFS_SKIP_SMUDGE=1 (skips large binary/lfs files)
# then pip-installing only the Python source.
#
# PyPI deps (transformers, datasets, huggingface_hub, etc.) are already installed
# by the previous 00_install_and_download.sh job.
#
# Model weights are downloaded separately by 00d_download_model_weights.sh.
#
# Submit: sbatch scripts/00c_install_geneformer_pkg.sh

set -euo pipefail

ENV_PATH="/work/pi_cfallini_uri_edu/riccardo-s_uri_edu-conda/.conda/envs/myenv"
PIP="${ENV_PATH}/bin/pip"
PYTHON="${ENV_PATH}/bin/python"
CLONE_DIR="/work/pi_cfallini_uri_edu/Geneformer_TBI/models/geneformer_src"

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
echo "Next: sbatch scripts/00d_download_model_weights.sh"
