#!/bin/bash
# Submit from the repository root. The log paths below are relative, because
# Slurm does not expand variables in #SBATCH lines, and the folder must exist:
#   mkdir -p results/02_geneformer_in_silico_screen/logs
#   sbatch 02_geneformer_in_silico_screen/00_environment/00_install_geneformer_env.sh
# Site-specific: partition, GPU constraint, node exclusion, array range, CPUs, memory and wall time are those of the original run (URI Unity cluster); adapt them to your scheduler.
#SBATCH --job-name=gf_install
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=results/02_geneformer_in_silico_screen/logs/install_geneformer_%j.out
#SBATCH --error=results/02_geneformer_in_silico_screen/logs/install_geneformer_%j.err

# Install Geneformer V2 and all perturbation dependencies into the environment
# given by GF_PYTHON. The environment it was run in (myenv) already had: Python 3.10, torch 2.7.1+cu118, anndata 0.11.4, scipy, pandas, h5py.
# This script adds: geneformer, transformers, datasets, huggingface_hub, einops, scanpy.

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

echo "=== Python version in myenv ==="
${PYTHON} --version

echo ""
echo "=== Installing Geneformer V2 and dependencies ==="
${PIP} install --upgrade \
    "geneformer" \
    "transformers>=4.36.0" \
    "datasets>=2.14.0" \
    "huggingface_hub>=0.20.0" \
    "einops>=0.6.0" \
    "accelerate>=0.25.0" \
    "pyarrow>=12.0.0" \
    "scanpy>=1.9.6"

echo ""
echo "=== Verifying installation ==="
${PYTHON} -c "
import geneformer
import torch
import transformers
import datasets
print('geneformer:', geneformer.__file__)
print('torch:', torch.__version__)
print('CUDA available:', torch.cuda.is_available())
print('transformers:', transformers.__version__)
print('datasets:', datasets.__version__)
print('All imports OK')
"

echo ""
echo "=== Done. Geneformer environment ready. ==="
