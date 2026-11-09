#!/bin/bash
#SBATCH --job-name=gf_install
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/install_geneformer_%j.out
#SBATCH --error=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/install_geneformer_%j.err

# Install Geneformer V2 and all perturbation dependencies into myenv.
# myenv already has: Python 3.10, torch 2.7.1+cu118, anndata 0.11.4, scipy, pandas, h5py.
# This script adds: geneformer, transformers, datasets, huggingface_hub, einops, scanpy.

set -euo pipefail

CONDA_BASE="/work/pi_cfallini_uri_edu/riccardo-s_uri_edu-conda"
ENV_PATH="${CONDA_BASE}/.conda/envs/myenv"
PIP="${ENV_PATH}/bin/pip"

echo "=== Python version in myenv ==="
${ENV_PATH}/bin/python --version

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
${ENV_PATH}/bin/python -c "
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
