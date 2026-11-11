#!/bin/bash
# Submit from the repository root. The log paths below are relative, because
# Slurm does not expand variables in #SBATCH lines, and the folder must exist:
#   mkdir -p results/02_geneformer_in_silico_screen/logs
#   sbatch 02_geneformer_in_silico_screen/00_environment/00e_fix_transformers_deps.sh
# Site-specific: partition, GPU constraint, node exclusion, array range, CPUs, memory and wall time are those of the original run (URI Unity cluster); adapt them to your scheduler.
#SBATCH --job-name=gf_deps_fix
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=results/02_geneformer_in_silico_screen/logs/gf_deps_fix_%j.out
#SBATCH --error=results/02_geneformer_in_silico_screen/logs/gf_deps_fix_%j.err
#SBATCH --mail-type=END,FAIL
# No --mail-user here: pass your own address on the command line (sbatch --mail-user=<address> ...).

# Geneformer 0.1.0 requires transformers==4.46.x.
# We installed 5.8.0 which removed SpecialTokensMixin -> ImportError.
# This script pins transformers to 4.46.* and installs remaining geneformer deps.

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

echo "Job: ${SLURM_JOB_ID}  |  Node: ${SLURMD_NODENAME}  |  $(date)"

echo ""
echo "--- Pinning transformers to 4.46.x (required by geneformer 0.1.0) ---"
${PIP} install "transformers==4.46.*" "tokenizers>=0.20,<0.22"

echo ""
echo "--- Installing remaining geneformer dependencies ---"
${PIP} install \
    "peft>=0.11.1" \
    "bitsandbytes>=0.45.5" \
    "optuna>=3.6" \
    "optuna-integration>=3.6" \
    "ray>=2.6" \
    "tdigest>=0.5.2" \
    "tensorboard>=2.15" \
    "loompy>=3.0" \
    "hyperopt>=0.2"

echo ""
echo "--- Verifying all key geneformer imports ---"
${PYTHON} -c "
import geneformer
from geneformer import TranscriptomeTokenizer
from geneformer import InSilicoPerturber
from geneformer import InSilicoPerturberStats
import torch
import transformers
print('geneformer:           ', geneformer.__file__)
print('TranscriptomeTokenizer OK')
print('InSilicoPerturber     OK')
print('InSilicoPerturberStats OK')
print('torch version:        ', torch.__version__)
print('transformers version: ', transformers.__version__)
print('CUDA (CPU node):      ', torch.cuda.is_available())
print()
print('All imports successful.')
"

echo ""
echo "Done at $(date)"
echo "Next: sbatch 02_geneformer_in_silico_screen/00_environment/00d_download_model_weights.sh"
