#!/bin/bash
# Submit from the repository root. The log paths below are relative, because
# Slurm does not expand variables in #SBATCH lines, and the folder must exist:
#   mkdir -p results/02_geneformer_in_silico_screen/logs
#   sbatch 02_geneformer_in_silico_screen/06_confounders_and_robustness/run_cross_direction_network_v2.sh
# Without Slurm the Python script can be run directly (it takes a few seconds);
# ../run_section.sh does that.
# Site-specific: partition, CPUs, memory and wall time are those of the original run (URI Unity cluster); adapt them to your scheduler.
#SBATCH --job-name=cross_dir_v2
#SBATCH -p cpu,cpu-preempt
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH --output=results/02_geneformer_in_silico_screen/logs/cross_direction_network_v2_%j.out
#SBATCH --error=results/02_geneformer_in_silico_screen/logs/cross_direction_network_v2_%j.err
#SBATCH --mail-type=END,FAIL
# No --mail-user here: pass your own address on the command line (sbatch --mail-user=<address> ...).

# v2 cross-direction network analysis — 4 datasets (adds RESCUE_ovr to D1/D2/D3).
# Hits STRING + Enrichr libraries on disk; outputs to results/cross_direction_network_v2/.
# In this repository the outputs go to
# <results>/02_geneformer_in_silico_screen/cross_direction_network_v2/ and the
# inputs are read below $GF_WORK/results/ (stage 07 and stage 09 outputs).

set -euo pipefail
# Repository and section configuration (gf_env.sh documents every GF_* variable).
# Under sbatch this script runs from a spool copy, so the repository is taken
# from the submit directory, or from NANOAI_REPO if that is exported.
GF_REPO="${NANOAI_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
if [ ! -f "${GF_REPO}/02_geneformer_in_silico_screen/gf_env.sh" ]; then
    echo "ERROR: repository not found at ${GF_REPO}; submit from the repository root or export NANOAI_REPO." >&2
    exit 1
fi
source "${GF_REPO}/02_geneformer_in_silico_screen/gf_env.sh"

PYTHON="${GF_PYTHON}"
SCRIPT="${GF_SECTION_DIR}/06_confounders_and_robustness/cross_direction_network_analysis.py"

T_START=$(date +%s)
echo "cross-direction network v2 | Job ${SLURM_JOB_ID} | Start $(date)"
${PYTHON} "${SCRIPT}"
echo ""
echo "Finished $(date) | elapsed $(( ($(date +%s) - T_START) / 60 )) min"
