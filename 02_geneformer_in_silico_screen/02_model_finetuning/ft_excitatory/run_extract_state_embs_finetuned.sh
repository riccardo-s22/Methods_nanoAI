#!/bin/bash
# Submit from the repository root. The log paths below are relative, because
# Slurm does not expand variables in #SBATCH lines, and the folder must exist:
#   mkdir -p results/02_geneformer_in_silico_screen/logs
#   sbatch 02_geneformer_in_silico_screen/02_model_finetuning/ft_excitatory/run_extract_state_embs_finetuned.sh
# Site-specific: partition, GPU constraint, node exclusion, array range, CPUs, memory and wall time are those of the original run (URI Unity cluster); adapt them to your scheduler.
#SBATCH --job-name=gf_tbi_ftstate
#SBATCH -p uri-gpu
#SBATCH --gpus=1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 00:30:00
#SBATCH --output=results/02_geneformer_in_silico_screen/logs/finetune_stateembs_%j.out
#SBATCH --error=results/02_geneformer_in_silico_screen/logs/finetune_stateembs_%j.err

# Extract TBI/control state embeddings using the fine-tuned classifier
# checkpoint. Single pass over 6,025 excitatory cells; ~3-5 min on A100.

set -euo pipefail
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TOKENIZERS_PARALLELISM=false

# Repository and section configuration (gf_env.sh documents every GF_* variable).
# Under sbatch this script runs from a spool copy, so the repository is taken
# from the submit directory, or from NANOAI_REPO if that is exported.
GF_REPO="${NANOAI_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}}"
if [ ! -f "${GF_REPO}/02_geneformer_in_silico_screen/gf_env.sh" ]; then
    echo "ERROR: repository not found at ${GF_REPO}; submit from the repository root or export NANOAI_REPO." >&2
    exit 1
fi
source "${GF_REPO}/02_geneformer_in_silico_screen/gf_env.sh"

PYTHON="${GF_PYTHON}"
SCRIPT="${GF_SECTION_DIR}/02_model_finetuning/ft_excitatory/13_extract_state_embs_finetuned.py"

echo "Job: ${SLURM_JOB_ID} on ${SLURMD_NODENAME} at $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || true

${PYTHON} "${SCRIPT}"

echo "Finished: $(date)"
