#!/bin/bash
#SBATCH --job-name=gf_tbi_ftstate
#SBATCH -p uri-gpu
#SBATCH --gpus=1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 00:30:00
#SBATCH --output=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/finetune_stateembs_%j.out
#SBATCH --error=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/finetune_stateembs_%j.err

# Extract TBI/control state embeddings using the fine-tuned classifier
# checkpoint. Single pass over 6,025 excitatory cells; ~3-5 min on A100.

set -euo pipefail
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TOKENIZERS_PARALLELISM=false

ENV_PATH="/work/pi_cfallini_uri_edu/riccardo-s_uri_edu-conda/.conda/envs/myenv"
PYTHON="${ENV_PATH}/bin/python"
SCRIPT="/work/pi_cfallini_uri_edu/Geneformer_TBI/scripts/13_extract_state_embs_finetuned.py"

echo "Job: ${SLURM_JOB_ID} on ${SLURMD_NODENAME} at $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || true

${PYTHON} "${SCRIPT}"

echo "Finished: $(date)"
