#!/bin/bash
# Submit from the repository root. The log paths below are relative, because
# Slurm does not expand variables in #SBATCH lines, and the folder must exist:
#   mkdir -p results/02_geneformer_in_silico_screen/logs
#   sbatch 02_geneformer_in_silico_screen/02_model_finetuning/ft_excitatory/run_finetune_tbi_classifier_excitatory.sh
# Site-specific: partition, GPU constraint, node exclusion, array range, CPUs, memory and wall time are those of the original run (URI Unity cluster); adapt them to your scheduler.
#SBATCH --job-name=gf_tbi_ftune
#SBATCH -p uri-gpu
#SBATCH --gpus=1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 06:00:00
# uri-gpu only: V2 104M fine-tune at batch=12 needs >11 GB VRAM; the broader
# gpu pool can land on 11 GB GTX 1080 Ti nodes and OOMs in backward pass.
#SBATCH --output=results/02_geneformer_in_silico_screen/logs/finetune_tbi_excit_%j.out
#SBATCH --error=results/02_geneformer_in_silico_screen/logs/finetune_tbi_excit_%j.err
#SBATCH --mail-type=END,FAIL
# No --mail-user here: pass your own address on the command line (sbatch --mail-user=<address> ...).

# Fine-tune Geneformer-V2-104M for TBI vs. control cell classification on
# excitatory neurons with a donor-level held-out split.
# Held-out donors: control 502, TBI 6, TBI 19.

set -euo pipefail
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TOKENIZERS_PARALLELISM=false
export WANDB_DISABLED=true

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
SCRIPT="${GF_SECTION_DIR}/02_model_finetuning/ft_excitatory/12_finetune_tbi_classifier_excitatory.py"

T_START=$(date +%s)

echo "============================================================"
echo "Geneformer TBI — excitatory TBI classifier fine-tuning"
echo "Job:    ${SLURM_JOB_ID}"
echo "Node:   ${SLURMD_NODENAME}"
echo "Start:  $(date)"
echo "------------------------------------------------------------"
nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv,noheader 2>/dev/null || echo "nvidia-smi: N/A"
echo "============================================================"

${PYTHON} "${SCRIPT}"

echo ""
echo "------------------------------------------------------------"
echo "Finished at: $(date)"
echo "Elapsed: $(( ($(date +%s) - T_START) / 60 )) min"
nvidia-smi --query-gpu=name,memory.used,memory.free --format=csv,noheader 2>/dev/null || true
