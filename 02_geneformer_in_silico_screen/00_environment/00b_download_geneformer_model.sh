#!/bin/bash
#SBATCH --job-name=gf_download
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/download_model_%j.out
#SBATCH --error=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/download_model_%j.err

# Download Geneformer V2 pre-trained model weights from HuggingFace.
# Downloads the 12-layer (12L) pre-trained backbone used for in silico perturbation.
# Model destination: /work/pi_cfallini_uri_edu/Geneformer_TBI/models/Geneformer/

set -euo pipefail

ENV_PATH="/work/pi_cfallini_uri_edu/riccardo-s_uri_edu-conda/.conda/envs/myenv"
PYTHON="${ENV_PATH}/bin/python"
MODEL_DIR="/work/pi_cfallini_uri_edu/Geneformer_TBI/models"

echo "=== Downloading Geneformer V2 pre-trained models ==="
echo "Destination: ${MODEL_DIR}"

${PYTHON} - <<'PYEOF'
import os
from pathlib import Path
from huggingface_hub import snapshot_download

model_dir = Path("/work/pi_cfallini_uri_edu/Geneformer_TBI/models")

# Download the full Geneformer model repository (includes 6L and 12L pre-trained weights
# plus fine-tuned models for various tasks, and the gene token dictionaries).
print("Downloading ctheodoris/Geneformer ...")
local_path = snapshot_download(
    repo_id="ctheodoris/Geneformer",
    local_dir=str(model_dir / "Geneformer"),
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
echo "Model weights are in: ${MODEL_DIR}/Geneformer/"
echo ""
echo "Next steps:"
echo "  1. Check which pre-trained model directory to use for ISP"
echo "  2. Run pilot perturbation: sbatch scripts/run_pilot_perturb.sh"
