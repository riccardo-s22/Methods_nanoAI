#!/bin/bash
#SBATCH --job-name=gf_weights
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/gf_weights_%j.out
#SBATCH --error=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/gf_weights_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=riccardo-s@uri.edu

# Download Geneformer V2 pre-trained model weights via huggingface_hub snapshot_download.
# This properly handles git-lfs files without requiring git-lfs on the cluster.
#
# Requires: huggingface_hub installed (done by 00_install_and_download.sh).
# Submit after 00c completes: sbatch scripts/00d_download_model_weights.sh

set -euo pipefail

ENV_PATH="/work/pi_cfallini_uri_edu/riccardo-s_uri_edu-conda/.conda/envs/myenv"
PYTHON="${ENV_PATH}/bin/python"
MODEL_DIR="/work/pi_cfallini_uri_edu/Geneformer_TBI/models/Geneformer"

echo "Job: ${SLURM_JOB_ID}  |  Node: ${SLURMD_NODENAME}  |  $(date)"
echo "Downloading to: ${MODEL_DIR}"

${PYTHON} - <<'PYEOF'
from pathlib import Path
from huggingface_hub import snapshot_download

model_dir = Path("/work/pi_cfallini_uri_edu/Geneformer_TBI/models/Geneformer")
model_dir.mkdir(parents=True, exist_ok=True)

print("Downloading ctheodoris/Geneformer model weights ...")
print("(This includes pre-trained 6L and 12L backbone weights; may take 20-60 min)")

local_path = snapshot_download(
    repo_id="ctheodoris/Geneformer",
    local_dir=str(model_dir),
    ignore_patterns=["*.msgpack", "flax_model*", "tf_model*", "*.h5"],
)
print(f"\nDownloaded to: {local_path}")

print("\n=== Top-level contents ===")
for p in sorted(Path(local_path).iterdir()):
    if p.is_dir():
        n_files = len(list(p.rglob("*")))
        size    = sum(f.stat().st_size for f in p.rglob("*") if f.is_file())
        print(f"  [dir]  {p.name}/  ({n_files} files, {size/1024**2:.0f} MB)")
    elif p.is_file():
        print(f"  [file] {p.name}  ({p.stat().st_size/1024**2:.1f} MB)")

# Find and report the pre-trained backbone directories
print("\n=== Pre-trained backbone models found ===")
for cfg in sorted(Path(local_path).rglob("config.json")):
    d = cfg.parent
    if d == Path(local_path):
        continue
    size = sum(f.stat().st_size for f in d.rglob("*") if f.is_file())
    print(f"  {d.relative_to(local_path)}  ({size/1024**2:.0f} MB)")
PYEOF

echo ""
echo "Download complete at $(date)"
echo ""
echo "Next steps:"
echo "  1. python scripts/01_count_neuronal_genes.py"
echo "  2. sbatch scripts/run_pilot_perturb.sh"
