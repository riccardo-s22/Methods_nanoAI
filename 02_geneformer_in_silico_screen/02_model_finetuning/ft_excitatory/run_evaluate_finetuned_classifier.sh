#!/bin/bash
#SBATCH --job-name=gf_tbi_eval
#SBATCH -p uri-gpu
#SBATCH --gpus=1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 00:30:00
#SBATCH --output=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/finetune_eval_%j.out
#SBATCH --error=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/finetune_eval_%j.err

# Eval-only run for the previously-trained TBI classifier (training succeeded
# in job 57711067; only the eval step had a column-name bug).

set -euo pipefail
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TOKENIZERS_PARALLELISM=false
export WANDB_DISABLED=true

ENV_PATH="/work/pi_cfallini_uri_edu/riccardo-s_uri_edu-conda/.conda/envs/myenv"
PYTHON="${ENV_PATH}/bin/python"
SCRIPT="/work/pi_cfallini_uri_edu/Geneformer_TBI/scripts/12b_evaluate_finetuned_classifier.py"

echo "Job: ${SLURM_JOB_ID} on ${SLURMD_NODENAME} at $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || true

${PYTHON} "${SCRIPT}"

echo "Finished: $(date)"
