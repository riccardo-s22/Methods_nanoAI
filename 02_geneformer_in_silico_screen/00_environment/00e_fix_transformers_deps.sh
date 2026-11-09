#!/bin/bash
#SBATCH --job-name=gf_deps_fix
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/gf_deps_fix_%j.out
#SBATCH --error=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/gf_deps_fix_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=riccardo-s@uri.edu

# Geneformer 0.1.0 requires transformers==4.46.x.
# We installed 5.8.0 which removed SpecialTokensMixin -> ImportError.
# This script pins transformers to 4.46.* and installs remaining geneformer deps.

set -euo pipefail

ENV_PATH="/work/pi_cfallini_uri_edu/riccardo-s_uri_edu-conda/.conda/envs/myenv"
PIP="${ENV_PATH}/bin/pip"
PYTHON="${ENV_PATH}/bin/python"

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
echo "Next: sbatch scripts/00d_download_model_weights.sh"
