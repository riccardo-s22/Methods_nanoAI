#!/bin/bash
#SBATCH --job-name=gf_setup
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/gf_setup_%j.out
#SBATCH --error=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/gf_setup_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=riccardo-s@uri.edu

# Single-job setup: installs Geneformer + all deps, then downloads model weights.
# Steps run in order so huggingface_hub is available before the download.
#
# Uses the full path to myenv (no module load needed):
#   Python 3.10, torch 2.7.1+cu118, anndata 0.11.4 already present.
#
# Submit: sbatch scripts/00_install_and_download.sh

set -euo pipefail

ENV_PATH="/work/pi_cfallini_uri_edu/riccardo-s_uri_edu-conda/.conda/envs/myenv"
PIP="${ENV_PATH}/bin/pip"
PYTHON="${ENV_PATH}/bin/python"
MODEL_DIR="/work/pi_cfallini_uri_edu/Geneformer_TBI/models/Geneformer"

echo "============================================================"
echo "Geneformer TBI — Environment setup"
echo "Job: ${SLURM_JOB_ID}  |  Node: ${SLURMD_NODENAME}"
echo "Start: $(date)"
echo "============================================================"

echo ""
echo "--- Step 1: Python version ---"
${PYTHON} --version
${PIP} --version

echo ""
echo "--- Step 2: Install PyPI dependencies (includes huggingface_hub) ---"
${PIP} install --upgrade \
    "transformers>=4.36.0" \
    "datasets>=2.14.0" \
    "huggingface_hub>=0.20.0" \
    "einops>=0.6.0" \
    "accelerate>=0.25.0" \
    "pyarrow>=12.0.0" \
    "scanpy>=1.9.6"

echo ""
echo "--- Step 3: Install Geneformer from HuggingFace ---"
# geneformer is not on PyPI; install from the HF model/package repo
${PIP} install "git+https://huggingface.co/ctheodoris/Geneformer"

echo ""
echo "--- Step 4: Verify imports ---"
${PYTHON} -c "
import geneformer; print('geneformer:', geneformer.__file__)
import torch;       print('torch:', torch.__version__)
import transformers; print('transformers:', transformers.__version__)
import datasets;    print('datasets:', datasets.__version__)
print('CUDA available (CPU node, expected False):', torch.cuda.is_available())
print('All imports OK')
"

echo ""
echo "--- Step 5: Download Geneformer V2 model weights ---"
echo "Destination: ${MODEL_DIR}"
mkdir -p "${MODEL_DIR}"

${PYTHON} - <<'PYEOF'
from pathlib import Path
from huggingface_hub import snapshot_download

model_dir = Path("/work/pi_cfallini_uri_edu/Geneformer_TBI/models/Geneformer")
print("Downloading ctheodoris/Geneformer (this may take 20–40 minutes) ...")

local_path = snapshot_download(
    repo_id="ctheodoris/Geneformer",
    local_dir=str(model_dir),
    # Skip Flax/JAX weights and msgpack to save space
    ignore_patterns=["*.msgpack", "flax_model*", "tf_model*"],
)
print(f"Downloaded to: {local_path}")

# Report what landed
print("\n=== Model directory top-level ===")
for p in sorted(Path(local_path).iterdir()):
    if p.is_dir():
        size = sum(f.stat().st_size for f in p.rglob("*") if f.is_file())
        print(f"  [dir]  {p.name}/  ({size/1024**2:.0f} MB)")
    else:
        print(f"  [file] {p.name}  ({p.stat().st_size/1024**2:.1f} MB)")
PYEOF

echo ""
echo "============================================================"
echo "Setup complete at: $(date)"
echo ""
echo "Next steps:"
echo "  1. python scripts/01_count_neuronal_genes.py   (re-run to intersect with GF vocab)"
echo "  2. sbatch scripts/run_pilot_perturb.sh"
echo "============================================================"
