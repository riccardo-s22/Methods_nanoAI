#!/bin/bash
# Submit from the repository root. The log paths below are relative, because
# Slurm does not expand variables in #SBATCH lines, and the folder must exist:
#   mkdir -p results/02_geneformer_in_silico_screen/logs
#   sbatch 02_geneformer_in_silico_screen/00_environment/00d_download_model_weights.sh
# Site-specific: partition, GPU constraint, node exclusion, array range, CPUs, memory and wall time are those of the original run (URI Unity cluster); adapt them to your scheduler.
#SBATCH --job-name=gf_weights
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=results/02_geneformer_in_silico_screen/logs/gf_weights_%j.out
#SBATCH --error=results/02_geneformer_in_silico_screen/logs/gf_weights_%j.err
#SBATCH --mail-type=END,FAIL
# No --mail-user here: pass your own address on the command line (sbatch --mail-user=<address> ...).

# Download Geneformer V2 pre-trained model weights via huggingface_hub snapshot_download.
# This properly handles git-lfs files without requiring git-lfs on the cluster.
#
# Requires: huggingface_hub installed (done by 00_install_and_download.sh).
# Submit after 00c and 00e complete.

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

echo "Job: ${SLURM_JOB_ID}  |  Node: ${SLURMD_NODENAME}  |  $(date)"
echo "Downloading to: ${MODEL_DIR}"

${PYTHON} - <<'PYEOF'
import os
from pathlib import Path
from huggingface_hub import snapshot_download

model_dir = Path(os.environ["GF_MODEL_DIR"])
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
