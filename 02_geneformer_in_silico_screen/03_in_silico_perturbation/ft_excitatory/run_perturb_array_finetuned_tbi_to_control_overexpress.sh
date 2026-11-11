#!/bin/bash
# Submit from the repository root. The log paths below are relative, because
# Slurm does not expand variables in #SBATCH lines, and the folder must exist:
#   mkdir -p results/02_geneformer_in_silico_screen/logs
#   sbatch 02_geneformer_in_silico_screen/03_in_silico_perturbation/ft_excitatory/run_perturb_array_finetuned_tbi_to_control_overexpress.sh
# Site-specific: partition, GPU constraint, node exclusion, array range, CPUs, memory and wall time are those of the original run (URI Unity cluster); adapt them to your scheduler.
#SBATCH --job-name=gf_ft_ovr_rescue
#SBATCH -p uri-gpu,gpu
#SBATCH --constraint=a100-80g
#SBATCH --gpus=1
#SBATCH --exclude=gpu023
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH --array=0-36
#SBATCH --output=results/02_geneformer_in_silico_screen/logs/perturb_ft_ovr_rescue_%A_%a.out
#SBATCH --error=results/02_geneformer_in_silico_screen/logs/perturb_ft_ovr_rescue_%A_%a.err
#SBATCH --mail-type=END,FAIL
# No --mail-user here: pass your own address on the command line (sbatch --mail-user=<address> ...).

# Fine-tuned TBI → control RESCUE perturbation, OVEREXPRESS mode — excitatory.
# Top hits = candidate therapeutic-rescue genes (gain-of-function rescue).
# Same chunking as delete (37 chunks for 2,047 TBI cells).

set -euo pipefail
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export PERTURB_TYPE=overexpress

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
SCRIPT="${GF_SECTION_DIR}/03_in_silico_perturbation/ft_excitatory/perturb_cells_finetuned_tbi_to_control.py"

T_START=$(date +%s)
echo "============================================================"
echo "Geneformer TBI — fine-tuned TBI→control RESCUE (OVEREXPRESS)"
echo "Array job: ${SLURM_ARRAY_JOB_ID}  task: ${SLURM_ARRAY_TASK_ID}"
echo "Node:      ${SLURMD_NODENAME}"
echo "Start:     $(date)"
echo "------------------------------------------------------------"
nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv,noheader 2>/dev/null || echo "nvidia-smi: N/A"
echo "============================================================"

${PYTHON} "${SCRIPT}"

echo "------------------------------------------------------------"
echo "Task ${SLURM_ARRAY_TASK_ID} finished at: $(date)"
echo "Elapsed: $(( ($(date +%s) - T_START) / 60 )) min"
