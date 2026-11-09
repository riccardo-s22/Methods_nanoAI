#!/bin/bash
#SBATCH --job-name=gf_tbi_ftune
#SBATCH -p uri-gpu
#SBATCH --gpus=1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 06:00:00
# uri-gpu only: V2 104M fine-tune at batch=12 needs >11 GB VRAM; the broader
# gpu pool can land on 11 GB GTX 1080 Ti nodes and OOMs in backward pass.
#SBATCH --output=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/finetune_tbi_excit_%j.out
#SBATCH --error=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/finetune_tbi_excit_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=riccardo-s@uri.edu

# Fine-tune Geneformer-V2-104M for TBI vs. control cell classification on
# excitatory neurons with a donor-level held-out split.
# Held-out donors: control 502, TBI 6, TBI 19.

set -euo pipefail
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TOKENIZERS_PARALLELISM=false
export WANDB_DISABLED=true

ENV_PATH="/work/pi_cfallini_uri_edu/riccardo-s_uri_edu-conda/.conda/envs/myenv"
PYTHON="${ENV_PATH}/bin/python"
SCRIPT="/work/pi_cfallini_uri_edu/Geneformer_TBI/scripts/12_finetune_tbi_classifier_excitatory.py"

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
