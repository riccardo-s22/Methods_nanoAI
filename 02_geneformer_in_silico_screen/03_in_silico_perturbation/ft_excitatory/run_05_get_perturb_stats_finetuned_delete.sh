#!/bin/bash
#SBATCH --job-name=ft_stats_del
#SBATCH -p cpu,cpu-preempt
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --output=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/05_stats_finetuned_delete_%j.out
#SBATCH --error=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/05_stats_finetuned_delete_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=riccardo-s@uri.edu

# Aggregate fine-tuned control->TBI delete perturbation pickles into per-gene
# stats. Reads 10,452 pickles from chunks_delete/ and produces:
#   stats_delete/excit_finetuned_delete.csv                          (ISP raw)
#   stats_delete/excit_finetuned_delete_unfiltered_all_gene_ranked.tsv
#   stats_delete/excit_finetuned_delete_filtered_high_confidence_gene_ranked.tsv

set -euo pipefail
ENV_PATH="/work/pi_cfallini_uri_edu/riccardo-s_uri_edu-conda/.conda/envs/myenv"
PYTHON="${ENV_PATH}/bin/python"
SCRIPT="/work/pi_cfallini_uri_edu/Geneformer_TBI/scripts/05_get_perturb_stats_finetuned.py"

T_START=$(date +%s)
echo "GARZA stage 4 — aggregate delete stats | Job ${SLURM_JOB_ID} | Start $(date)"
${PYTHON} "${SCRIPT}" --mode delete
echo "Finished $(date) | elapsed $(( ($(date +%s) - T_START) / 60 )) min"
