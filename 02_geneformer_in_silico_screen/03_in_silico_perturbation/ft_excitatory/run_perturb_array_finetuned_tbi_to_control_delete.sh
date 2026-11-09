#!/bin/bash
#SBATCH --job-name=gf_ft_del_rescue
#SBATCH -p uri-gpu,gpu
#SBATCH --constraint=a100-80g
#SBATCH --gpus=1
#SBATCH --exclude=gpu023
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH --array=0-36
#SBATCH --output=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/perturb_ft_del_rescue_%A_%a.out
#SBATCH --error=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/perturb_ft_del_rescue_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=riccardo-s@uri.edu

# Fine-tuned TBI → control RESCUE perturbation, DELETE mode — excitatory neurons.
# Perturbs 2,047 TBI cells; measures cosine shift toward control centroid.
# Top hits = candidate therapeutic-rescue genes (loss-of-function rescue).
# 2,047 cells / 56 per chunk = 37 chunks (array 0-36; chunk 36 has 31 cells).
# Same A100-80g constraint as control→TBI screen (~5.7 min/cell baseline).

set -euo pipefail
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export PERTURB_TYPE=delete

ENV_PATH="/work/pi_cfallini_uri_edu/riccardo-s_uri_edu-conda/.conda/envs/myenv"
PYTHON="${ENV_PATH}/bin/python"
SCRIPT="/work/pi_cfallini_uri_edu/Geneformer_TBI/scripts/perturb_cells_finetuned_tbi_to_control.py"

T_START=$(date +%s)
echo "============================================================"
echo "Geneformer TBI — fine-tuned TBI→control RESCUE (DELETE)"
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
