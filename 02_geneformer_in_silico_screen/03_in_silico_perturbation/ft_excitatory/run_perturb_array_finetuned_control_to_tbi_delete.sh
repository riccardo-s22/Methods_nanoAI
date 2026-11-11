#!/bin/bash
# Submit from the repository root. The log paths below are relative, because
# Slurm does not expand variables in #SBATCH lines, and the folder must exist:
#   mkdir -p results/02_geneformer_in_silico_screen/logs
#   sbatch 02_geneformer_in_silico_screen/03_in_silico_perturbation/ft_excitatory/run_perturb_array_finetuned_control_to_tbi_delete.sh
# Site-specific: partition, GPU constraint, node exclusion, array range, CPUs, memory and wall time are those of the original run (URI Unity cluster); adapt them to your scheduler.
#SBATCH --job-name=gf_ft_del
#SBATCH -p uri-gpu,gpu
#SBATCH --constraint=bf16
#SBATCH --gpus=1
#SBATCH --exclude=gpu023
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH --array=0-71
#SBATCH --output=results/02_geneformer_in_silico_screen/logs/perturb_ft_del_%A_%a.out
#SBATCH --error=results/02_geneformer_in_silico_screen/logs/perturb_ft_del_%A_%a.err
#SBATCH --mail-type=END,FAIL
# No --mail-user here: pass your own address on the command line (sbatch --mail-user=<address> ...).

# Fine-tuned control→TBI perturbation, DELETE mode — excitatory neurons.
# Array task N processes control cells [N*56 .. N*56+55] of the control-filtered
# excitatory dataset (3,978 cells total → 72 chunks, array 0-71; chunk 71 has 2 cells).
#
# Empirical baseline from the zero-shot pipeline: ~434 s/cell on A100-40GB.
# Expected per-task: 56 cells × 434 s + ~7 min startup ≈ 6.9 h → 12 h wall.
#
# --constraint=bf16 restricts to Ampere+ GPUs (A100/L40S/H100). The original
# submission (57730381) let tasks land on gypsum-gpu* 1080 Ti / 2080 Ti nodes,
# which run ~2–3× slower and time out at the 12 h wall with no pickles emitted
# (28/72 tasks failed this way). bf16 is the cleanest Slurm feature filter
# spanning both partitions.
#
# Prerequisites:
#   - tokenized_paperQC_nGene1000/excitatory_neuron/*.dataset                                ✓
#   - results/finetune_tbi_excitatory/trained_model_path.txt → fine-tuned checkpoint        ✓
#   - results/perturb_neurons_finetuned_control_to_tbi/state_embs/excitatory_finetuned_state_embs.pkl  ✓
#
# After all 72 chunks complete:
#   sbatch 02_geneformer_in_silico_screen/03_in_silico_perturbation/ft_excitatory/run_05_get_perturb_stats_finetuned_delete.sh

set -euo pipefail
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export PERTURB_TYPE=delete

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
SCRIPT="${GF_SECTION_DIR}/03_in_silico_perturbation/ft_excitatory/perturb_cells_finetuned_control_to_tbi.py"

T_START=$(date +%s)

echo "============================================================"
echo "Geneformer TBI — fine-tuned control→TBI perturbation (DELETE)"
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
echo "Output:  ${GF_WORK}/results/perturb_neurons_finetuned_control_to_tbi/chunks_delete/"
nvidia-smi --query-gpu=name,memory.used,memory.free --format=csv,noheader 2>/dev/null || true
