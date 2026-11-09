#!/bin/bash
#SBATCH --job-name=gf_ft_oe
#SBATCH -p uri-gpu,gpu
#SBATCH --constraint=bf16
#SBATCH --gpus=1
#SBATCH --exclude=gpu023
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH --array=0-71
#SBATCH --output=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/perturb_ft_oe_%A_%a.out
#SBATCH --error=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/perturb_ft_oe_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=riccardo-s@uri.edu

# Fine-tuned control→TBI perturbation, OVEREXPRESS mode — excitatory neurons.
# Array task N processes control cells [N*56 .. N*56+55] of the control-filtered
# excitatory dataset (3,978 cells total → 72 chunks, array 0-71; chunk 71 has 2 cells).
#
# Empirical baseline from the zero-shot pipeline: ~434 s/cell on A100-40GB.
# Expected per-task: 56 cells × 434 s + ~7 min startup ≈ 6.9 h → 12 h wall.
#
# Overexpress mode moves each gene to the front of the rank-value encoding rather
# than deleting it from the encoding — same compute envelope as delete mode.
# (Geneformer's perturb_type literal is "overexpress"; this is the gene-activation
# modality from Theodoris 2023.)
#
# --constraint=bf16 restricts to Ampere+ GPUs (A100/L40S/H100). The original
# submission (57738575) let tasks land on gypsum-gpu* 1080 Ti / 2080 Ti nodes,
# which run ~2–3× slower and time out at the 12 h wall with no pickles emitted
# (25 confirmed TIMEOUTs + 10 doomed RUNNING tasks at last check). bf16 is the
# cleanest Slurm feature filter spanning both partitions.
#
# Prerequisites:
#   - tokenized_paperQC_nGene1000/excitatory_neuron/*.dataset                                ✓
#   - results/finetune_tbi_excitatory/trained_model_path.txt → fine-tuned checkpoint        ✓
#   - results/perturb_neurons_finetuned_control_to_tbi/state_embs/excitatory_finetuned_state_embs.pkl  ✓
#
# Submit:
#   sbatch scripts/run_perturb_array_finetuned_control_to_tbi_overexpress.sh
#
# After all 72 chunks complete:
#   python scripts/05_get_perturb_stats_finetuned.py --mode overexpress

set -euo pipefail
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export PERTURB_TYPE=overexpress

ENV_PATH="/work/pi_cfallini_uri_edu/riccardo-s_uri_edu-conda/.conda/envs/myenv"
PYTHON="${ENV_PATH}/bin/python"
SCRIPT="/work/pi_cfallini_uri_edu/Geneformer_TBI/scripts/perturb_cells_finetuned_control_to_tbi.py"

T_START=$(date +%s)

echo "============================================================"
echo "Geneformer TBI — fine-tuned control→TBI perturbation (OVEREXPRESS)"
echo "Array job: ${SLURM_ARRAY_JOB_ID}  task: ${SLURM_ARRAY_TASK_ID}"
echo "Node:      ${SLURMD_NODENAME}"
echo "Start:     $(date)"
echo "------------------------------------------------------------"
nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv,noheader 2>/dev/null || echo "nvidia-smi: N/A"
echo "============================================================"

${PYTHON} -c "
import torch, geneformer
print('torch:', torch.__version__, '| CUDA:', torch.cuda.is_available())
if torch.cuda.is_available():
    p = torch.cuda.get_device_properties(0)
    print('GPU:', p.name, '| VRAM:', round(p.total_memory/1e9,1), 'GB')
"

echo ""
${PYTHON} "${SCRIPT}"

echo ""
echo "------------------------------------------------------------"
echo "Task ${SLURM_ARRAY_TASK_ID} finished at: $(date)"
echo "Elapsed: $(( ($(date +%s) - T_START) / 60 )) min"
echo "Output:  /work/pi_cfallini_uri_edu/Geneformer_TBI/results/perturb_neurons_finetuned_control_to_tbi/chunks_overexpress/"
nvidia-smi --query-gpu=name,memory.used,memory.free --format=csv,noheader 2>/dev/null || true
