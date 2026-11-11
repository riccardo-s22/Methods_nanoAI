#!/bin/bash
# Submit from the repository root. The log paths below are relative, because
# Slurm does not expand variables in #SBATCH lines, and the folder must exist:
#   mkdir -p results/02_geneformer_in_silico_screen/logs
#   sbatch 02_geneformer_in_silico_screen/03_in_silico_perturbation/ft_excitatory/run_05_get_perturb_stats_finetuned_overexpress.sh
# Site-specific: partition, GPU constraint, node exclusion, array range, CPUs, memory and wall time are those of the original run (URI Unity cluster); adapt them to your scheduler.
#SBATCH --job-name=ft_stats_over
#SBATCH -p cpu,cpu-preempt
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --output=results/02_geneformer_in_silico_screen/logs/05_stats_finetuned_overexpress_%j.out
#SBATCH --error=results/02_geneformer_in_silico_screen/logs/05_stats_finetuned_overexpress_%j.err
#SBATCH --mail-type=END,FAIL
# No --mail-user here: pass your own address on the command line (sbatch --mail-user=<address> ...).

# Aggregate fine-tuned control->TBI OVEREXPRESS perturbation pickles into
# per-gene stats. Mirrors the delete launcher exactly; reads 10,452 pickles
# from chunks_overexpress/ and produces stats_overexpress/*.

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
SCRIPT="${GF_SECTION_DIR}/03_in_silico_perturbation/ft_excitatory/05_get_perturb_stats_finetuned.py"

T_START=$(date +%s)
echo "GARZA stage 4 — aggregate overexpress stats | Job ${SLURM_JOB_ID} | Start $(date)"
${PYTHON} "${SCRIPT}" --mode overexpress
echo "Finished $(date) | elapsed $(( ($(date +%s) - T_START) / 60 )) min"
