#!/bin/bash
#SBATCH --job-name=ft_stats_over
#SBATCH -p cpu,cpu-preempt
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --output=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/05_stats_finetuned_overexpress_%j.out
#SBATCH --error=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/05_stats_finetuned_overexpress_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=riccardo-s@uri.edu

# Aggregate fine-tuned control->TBI OVEREXPRESS perturbation pickles into
# per-gene stats. Mirrors the delete launcher exactly; reads 10,452 pickles
# from chunks_overexpress/ and produces stats_overexpress/*.

set -euo pipefail
ENV_PATH="/work/pi_cfallini_uri_edu/riccardo-s_uri_edu-conda/.conda/envs/myenv"
PYTHON="${ENV_PATH}/bin/python"
SCRIPT="/work/pi_cfallini_uri_edu/Geneformer_TBI/scripts/05_get_perturb_stats_finetuned.py"

T_START=$(date +%s)
echo "GARZA stage 4 — aggregate overexpress stats | Job ${SLURM_JOB_ID} | Start $(date)"
${PYTHON} "${SCRIPT}" --mode overexpress
echo "Finished $(date) | elapsed $(( ($(date +%s) - T_START) / 60 )) min"
