#!/bin/bash
#SBATCH --job-name=gf_ft_ovr_rescue
#SBATCH -p uri-gpu,gpu
#SBATCH --constraint=a100-80g
#SBATCH --gpus=1
#SBATCH --exclude=gpu023
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH --array=0-36
#SBATCH --output=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/perturb_ft_ovr_rescue_%A_%a.out
#SBATCH --error=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/perturb_ft_ovr_rescue_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=riccardo-s@uri.edu

# Fine-tuned TBI → control RESCUE perturbation, OVEREXPRESS mode — excitatory.
# Top hits = candidate therapeutic-rescue genes (gain-of-function rescue).
# Same chunking as delete (37 chunks for 2,047 TBI cells).

set -euo pipefail
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export PERTURB_TYPE=overexpress

ENV_PATH="/work/pi_cfallini_uri_edu/riccardo-s_uri_edu-conda/.conda/envs/myenv"
PYTHON="${ENV_PATH}/bin/python"
SCRIPT="/work/pi_cfallini_uri_edu/Geneformer_TBI/scripts/perturb_cells_finetuned_tbi_to_control.py"

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
