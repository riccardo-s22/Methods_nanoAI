#!/bin/bash
#SBATCH --job-name=ft_pipe_rsc_ovr
#SBATCH -p cpu,cpu-preempt
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --output=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/05_07_09_rescue_overexpress_%j.out
#SBATCH --error=/work/pi_cfallini_uri_edu/Geneformer_TBI/logs/05_07_09_rescue_overexpress_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=riccardo-s@uri.edu

# RESCUE (TBI->control) OVEREXPRESS pipeline — fine-tuned excitatory model.
# Chains stages 05 (per-gene stats), 07 (biological candidate filter),
# 09 (pathway enrichment) against chunks_overexpress/ (all 37 chunks complete).

set -euo pipefail
ENV_PATH="/work/pi_cfallini_uri_edu/riccardo-s_uri_edu-conda/.conda/envs/myenv"
PYTHON="${ENV_PATH}/bin/python"
SCRIPTS="/work/pi_cfallini_uri_edu/Geneformer_TBI/scripts"
CELLTYPE="excit_rescue"
MODE="overexpress"

cd /work/pi_cfallini_uri_edu/Geneformer_TBI

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
${PYTHON} "${SCRIPTS}/09_pathway_enrichment_finetuned.py" --mode "${MODE}" --celltype "${CELLTYPE}"

echo ""
echo "Finished $(date) | elapsed $(( ($(date +%s) - T_START) / 60 )) min"
