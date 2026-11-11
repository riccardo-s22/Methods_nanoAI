#!/bin/bash
# Submit from the repository root. The log paths below are relative, because
# Slurm does not expand variables in #SBATCH lines, and the folder must exist:
#   mkdir -p results/02_geneformer_in_silico_screen/logs
#   sbatch 02_geneformer_in_silico_screen/03_in_silico_perturbation/ft_excitatory/run_05_07_09_finetuned_rescue_overexpress.sh
# Site-specific: partition, GPU constraint, node exclusion, array range, CPUs, memory and wall time are those of the original run (URI Unity cluster); adapt them to your scheduler.
#SBATCH --job-name=ft_pipe_rsc_ovr
#SBATCH -p cpu,cpu-preempt
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --output=results/02_geneformer_in_silico_screen/logs/05_07_09_rescue_overexpress_%j.out
#SBATCH --error=results/02_geneformer_in_silico_screen/logs/05_07_09_rescue_overexpress_%j.err
#SBATCH --mail-type=END,FAIL
# No --mail-user here: pass your own address on the command line (sbatch --mail-user=<address> ...).

# RESCUE (TBI->control) OVEREXPRESS pipeline — fine-tuned excitatory model.
# Chains stages 05 (per-gene stats), 07 (biological candidate filter),
# 09 (pathway enrichment) against chunks_overexpress/ (all 37 chunks complete).

set -euo pipefail
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
# The three stages sat in one flat scripts/ folder; here 05 and 07 are in this
# folder and 09 is in 04_pathway_enrichment.
SCRIPTS="${GF_SECTION_DIR}/03_in_silico_perturbation/ft_excitatory"
SCRIPTS_09="${GF_SECTION_DIR}/04_pathway_enrichment"
CELLTYPE="excit_rescue"
MODE="overexpress"

cd "${GF_WORK}"

T_START=$(date +%s)
echo "RESCUE ovr pipeline | Job ${SLURM_JOB_ID} | Start $(date)"
echo "celltype=${CELLTYPE}  mode=${MODE}"
echo "=============================================="
echo "Stage 05: per-gene stats"
echo "=============================================="
${PYTHON} "${SCRIPTS}/05_get_perturb_stats_finetuned.py" --mode "${MODE}" --celltype "${CELLTYPE}"

echo ""
echo "=============================================="
echo "Stage 07: biological candidate filter"
echo "=============================================="
${PYTHON} "${SCRIPTS}/07_biological_candidate_filter_finetuned.py" --mode "${MODE}" --celltype "${CELLTYPE}"

echo ""
echo "=============================================="
echo "Stage 09: pathway enrichment"
echo "=============================================="
${PYTHON} "${SCRIPTS_09}/09_pathway_enrichment_finetuned.py" --mode "${MODE}" --celltype "${CELLTYPE}"

echo ""
echo "Finished $(date) | elapsed $(( ($(date +%s) - T_START) / 60 )) min"
