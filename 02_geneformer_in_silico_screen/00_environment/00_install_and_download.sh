#!/bin/bash
# Submit from the repository root. The log paths below are relative, because
# Slurm does not expand variables in #SBATCH lines, and the folder must exist:
#   mkdir -p results/02_geneformer_in_silico_screen/logs
#   sbatch 02_geneformer_in_silico_screen/00_environment/00_install_and_download.sh
# Site-specific: partition, GPU constraint, node exclusion, array range, CPUs, memory and wall time are those of the original run (URI Unity cluster); adapt them to your scheduler.
#SBATCH --job-name=gf_setup
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=results/02_geneformer_in_silico_screen/logs/gf_setup_%j.out
#SBATCH --error=results/02_geneformer_in_silico_screen/logs/gf_setup_%j.err
#SBATCH --mail-type=END,FAIL
# No --mail-user here: pass your own address on the command line (sbatch --mail-user=<address> ...).

# Single-job setup: installs Geneformer + all deps, then downloads model weights.
# Steps run in order so huggingface_hub is available before the download.
#
# Installs into the interpreter given by GF_PYTHON (no module load needed).
# The environment it was run in already held:
#   Python 3.10, torch 2.7.1+cu118, anndata 0.11.4.

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
MODEL_DIR="${GF_MODEL_DIR}"

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
import os
from pathlib import Path
from huggingface_hub import snapshot_download

model_dir = Path(os.environ["GF_MODEL_DIR"])
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
