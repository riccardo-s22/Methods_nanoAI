#!/bin/bash
#SBATCH --job-name=cross_dir_v2
#SBATCH -p cpu,cpu-preempt
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH --output=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/cross_direction_network_v2_%j.out
#SBATCH --error=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/cross_direction_network_v2_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=riccardo-s@uri.edu

# v2 cross-direction network analysis — 4 datasets (adds RESCUE_ovr to D1/D2/D3).
# Hits STRING + Enrichr libraries on disk; outputs to results/cross_direction_network_v2/.

set -euo pipefail
ENV_PATH="/work/pi_cfallini_uri_edu/riccardo-s_uri_edu-conda/.conda/envs/myenv"
PYTHON="${ENV_PATH}/bin/python"
SCRIPT="/work/pi_cfallini_uri_edu/Geneformer_TBI/scripts/cross_direction_network_analysis.py"

cd /work/pi_cfallini_uri_edu/Geneformer_TBI

T_START=$(date +%s)
echo "cross-direction network v2 | Job ${SLURM_JOB_ID} | Start $(date)"
${PYTHON} "${SCRIPT}"
echo ""
echo "Finished $(date) | elapsed $(( ($(date +%s) - T_START) / 60 )) min"
