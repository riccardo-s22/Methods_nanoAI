#!/bin/bash
# Submit from the repository root. The log paths below are relative, because
# Slurm does not expand variables in #SBATCH lines, and the folder must exist:
#   mkdir -p results/02_geneformer_in_silico_screen/logs
#   sbatch 02_geneformer_in_silico_screen/00_environment/00b_download_geneformer_model.sh
# Site-specific: partition, GPU constraint, node exclusion, array range, CPUs, memory and wall time are those of the original run (URI Unity cluster); adapt them to your scheduler.
#SBATCH --job-name=gf_download
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=results/02_geneformer_in_silico_screen/logs/download_model_%j.out
#SBATCH --error=results/02_geneformer_in_silico_screen/logs/download_model_%j.err

# Download Geneformer V2 pre-trained model weights from HuggingFace.
# Downloads the 12-layer (12L) pre-trained backbone used for in silico perturbation.
# Model destination: $GF_MODEL_DIR (default $GF_WORK/models/Geneformer)

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
MODEL_DIR="${GF_MODEL_DIR}"

echo "=== Downloading Geneformer V2 pre-trained models ==="
echo "Destination: ${MODEL_DIR}"

${PYTHON} - <<'PYEOF'
import os
from pathlib import Path
from huggingface_hub import snapshot_download

model_dir = Path(os.environ["GF_MODEL_DIR"])

# Download the full Geneformer model repository (includes 6L and 12L pre-trained weights
# plus fine-tuned models for various tasks, and the gene token dictionaries).
print("Downloading ctheodoris/Geneformer ...")
local_path = snapshot_download(
    repo_id="ctheodoris/Geneformer",
    local_dir=str(model_dir),
    ignore_patterns=["*.msgpack", "*.h5", "flax_model*"],
)
print(f"Downloaded to: {local_path}")

# List what we got
print("\n=== Model directory contents ===")
for p in sorted(Path(local_path).rglob("*"))[:80]:
    if p.is_file():
        size_mb = p.stat().st_size / 1024**2
        print(f"  {p.relative_to(local_path)}  ({size_mb:.1f} MB)")
PYEOF

echo ""
echo "=== Download complete ==="
echo "Model weights are in: ${MODEL_DIR}/"
echo ""
echo "Next steps:"
echo "  1. Check which pre-trained model directory to use for ISP"
echo "  2. Run pilot perturbation: sbatch scripts/run_pilot_perturb.sh"
