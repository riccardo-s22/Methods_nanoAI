#!/usr/bin/env bash
# Results 2 - Geneformer in-silico perturbation screen (cross-model consensus panel, Supp. Figs. 5-6).
# Re-runs every step of this section that can run on an ordinary CPU machine from the tables
# shipped in data/, in the order of README.md, and prints one SKIP line for every step that cannot.
#
#     bash 02_geneformer_in_silico_screen/run_section.sh [--offline] [--with-project-tree] [--only REGEX]
#     bash 02_geneformer_in_silico_screen/run_section.sh --print-cluster-plan
#
# Works from any working directory. Nothing is read from or written to the current directory:
#   inputs    <repo>/data/02_geneformer_in_silico_screen      (NANOAI_DATA overrides <repo>/data)
#   outputs   <repo>/results/02_geneformer_in_silico_screen   (NANOAI_RESULTS overrides <repo>/results)
#             and, for stages 07 and 09 only, the pipeline working tree $GF_WORK/results/
#             (default <repo>/work/Geneformer_TBI, git-ignored), where those stages have always
#             written and where the cross-direction analysis reads them
#   logs      <repo>/results/02_geneformer_in_silico_screen/logs/<step>.log
#
#   --offline             skip the five steps that query a web service (stage 09: Enrichr;
#                         cross-direction analysis: STRING). Also: GF_OFFLINE=1. By default they run.
#   --with-project-tree   also run the steps that need files which are not shipped in data/, when
#                         those files are found below $GF_WORK / $GF_JOINTPG / $GF_MECHANO (a copy
#                         of the original project trees, or the output of the cluster stages).
#                         Each such step still prints SKIP, naming the first missing file, when
#                         its inputs are absent. Also: GF_WITH_PROJECT_TREE=1.
#   --only REGEX          run only the steps whose name matches the extended regular expression,
#                         e.g. --only 'heatmap'  or  --only '^(seed|07_)'. Step names are the names
#                         in the PASS/FAIL/SKIP lines and of the log files.
#   --print-cluster-plan  print the upstream (cluster / GPU) half of the section as commands in
#                         run order, without executing anything, and exit.
#
# Interpreters: the default steps use ${PYTHON} (shared/env.sh, default python3) and need numpy,
# pandas, scipy, matplotlib, requests and gseapy. The --with-project-tree steps use ${GF_PREP_PYTHON}
# (anndata, scanpy, statsmodels) and ${GF_PYTHON} (geneformer, pyarrow, umap-learn); see gf_env.sh.
#
# Every step prints one PASS / FAIL / SKIP line; all of them are repeated in a summary at the end.
# The exit status is non-zero if any step failed. A failing step does not stop the run.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../shared/env.sh
source "${HERE}/../shared/env.sh"
# Section configuration: GF_WORK, GF_JOINTPG, GF_MECHANO, GF_PYTHON, GF_PREP_PYTHON, GF_DATA, GF_RESULTS.
# shellcheck source=gf_env.sh
source "${HERE}/gf_env.sh"

SECTION="02_geneformer_in_silico_screen"
OUT="${GF_RESULTS}"
LOGS="${OUT}/logs"
OFFLINE="${GF_OFFLINE:-0}"
WITH_TREE="${GF_WITH_PROJECT_TREE:-0}"
ONLY=""
PLAN=0

while [ $# -gt 0 ]; do
  case "$1" in
    --offline)             OFFLINE=1 ;;
    --with-project-tree)   WITH_TREE=1 ;;
    --print-cluster-plan)  PLAN=1 ;;
    --only)      [ $# -ge 2 ] || { echo "--only needs a regular expression" >&2; exit 2; }
                 ONLY="$2"; shift ;;
    --only=*)    ONLY="${1#--only=}" ;;
    -h|--help)   sed -n '2,34p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1 (use --help)" >&2; exit 2 ;;
  esac
  shift
done

# ------------------------------------------------------------------------------ the cluster plan
print_cluster_plan() {
  local S="${SECTION}"
  cat <<EOF
# Upstream half of ${S}: environment, dataset preparation, fine-tuning,
# in-silico perturbation. Nothing below is executed by run_section.sh. The commands are the
# original run order; partitions, constraints and limits in the Slurm headers are those of
# the original cluster and need adapting. Run everything from the repository root.

# ---- configuration (once per shell; every variable is documented in ${S}/gf_env.sh)
cd ${NANOAI_REPO}
export GF_CONDA_ENV=/path/to/conda/env            # Python 3.10 + CUDA PyTorch (original: Python 3.10.18, torch 2.7.1+cu118)
export GF_PREP_PYTHON=/path/to/prep/env/bin/python  # scanpy 1.11.5, harmonypy 2.0.0, scrublet 0.2.3, leidenalg, igraph
export GF_WORK=/path/to/Geneformer_TBI            # optional; default ${NANOAI_REPO}/work/Geneformer_TBI
source ${S}/gf_env.sh
mkdir -p results/${S}/logs               # Slurm writes its logs here
# sbatch forwards the exported variables. Add --mail-user=<address> to sbatch if wanted.
# No script creates the conda environment(s). envs/geneformer.yml is the recipe: one environment
# for both roles (then GF_PREP_PYTHON can stay unset), pins transcribed from the paper and not
# tested here; read its header first. The stage 00 jobs below install into that environment.

# ---- stage 00: environment (CPU partition, 4 CPUs, 16-32 GB, 1-3 h limits; network needed)
sbatch ${S}/00_environment/00_install_and_download.sh       # PyPI dependencies
sbatch ${S}/00_environment/00c_install_geneformer_pkg.sh    # clone ctheodoris/Geneformer without LFS, pip install --no-deps
sbatch ${S}/00_environment/00e_fix_transformers_deps.sh     # pins transformers 4.46.*
sbatch ${S}/00_environment/00d_download_model_weights.sh    # snapshot to \$GF_MODEL_DIR (2.8 GB)
# 00_install_geneformer_env.sh and 00b_download_geneformer_model.sh failed in the original run and are superseded.

# ---- stage 01: dataset preparation (CPU; about 16 GB RAM for the QC/annotation steps)
bash ${S}/01_dataset_preparation/fetch_geo_GSE209552.sh                           # 0.9 GB download; writes \$GF_GEO and reports/snrna_sample_manifest.txt
\$GF_PREP_PYTHON ${S}/01_dataset_preparation/parse_geo_series_metadata.py          # -> metadata/gse209552_sample_metadata_curated.tsv (see README: region caveat)
\$GF_PREP_PYTHON ${S}/01_dataset_preparation/build_geneformer_h5ad_from_10x.py     # -> h5ad/Garza_GSE209552_all_raw_counts_geneformer_ready.h5ad (24,621 nuclei)
\$GF_PREP_PYTHON ${S}/01_dataset_preparation/rebuild_with_paper_qc.py              # -> ..._paperQC_annotated.h5ad (first-pass QC, Scrublet, Harmony, Leiden, marker scoring)
\$GF_PREP_PYTHON ${S}/01_dataset_preparation/qc_sweep.py                           # -> reports/qc_sweep.tsv
\$GF_PREP_PYTHON ${S}/01_dataset_preparation/build_paper_matched_primary.py        # -> ..._paperQC_nGene1000_annotated.h5ad (15,172 nuclei) and h5ad/celltype_subsets_paperQC_nGene1000/
# manual step, not scripted: copy each subset into the folder the tokenizer reads
for pair in excitatory_neuron:excitatory_neurons inhibitory_neuron:inhibitory_neurons oligodendrocyte:oligodendrocyte OPC:OPC; do
  mkdir -p "\$GF_WORK/h5ad/geneformer_inputs_paperQC_nGene1000/\${pair%%:*}"
  cp "\$GF_WORK/h5ad/celltype_subsets_paperQC_nGene1000/\${pair##*:}_raw_geneformer_ready.h5ad" "\$GF_WORK/h5ad/geneformer_inputs_paperQC_nGene1000/\${pair%%:*}/"
done
\$GF_PYTHON ${S}/01_dataset_preparation/tokenize_paperQC_nGene1000.py              # -> tokenized_paperQC_nGene1000/<type>/*.dataset
# Superseded first attempts, kept as a record, not part of this sequence: organize_10x_from_geo_tar.sh,
# annotate_cell_types_de_novo.py, create_geneformer_subsets.py, rebuild_subsets_and_tokenize.py,
# validate_geneformer_h5ad.py, 01_count_neuronal_genes.py, apply_ngene1000_filter.py, refine_qc_to_paper.py.
# The last two write the same output names as build_paper_matched_primary.py: do not run them after it.

# ---- stage 02: FT fine-tuning (1 GPU with more than 11 GB memory)
sbatch ${S}/02_model_finetuning/ft_excitatory/run_finetune_tbi_classifier_excitatory.sh   # 8 CPUs, 64 GB, 6 h limit; original: 306 steps in 995 s
# only to evaluate an existing checkpoint (set GF_FT_MODEL to it): 4 CPUs, 32 GB, 30 min
sbatch ${S}/02_model_finetuning/ft_excitatory/run_evaluate_finetuned_classifier.sh
sbatch ${S}/02_model_finetuning/ft_excitatory/run_extract_state_embs_finetuned.sh         # 4 CPUs, 32 GB, 30 min -> state-embedding pickle

# ---- stage 03: in-silico perturbation (GPU arrays: 1 GPU, 8 CPUs, 32 GB, 12 h limit per task)
sbatch ${S}/03_in_silico_perturbation/ft_excitatory/run_perturb_array_finetuned_tbi_to_control_delete.sh        # rescue / delete,      --array=0-36
sbatch ${S}/03_in_silico_perturbation/ft_excitatory/run_perturb_array_finetuned_tbi_to_control_overexpress.sh   # rescue / overexpress, --array=0-36
sbatch ${S}/03_in_silico_perturbation/ft_excitatory/run_perturb_array_finetuned_control_to_tbi_delete.sh        # induction / delete,      --array=0-71
sbatch ${S}/03_in_silico_perturbation/ft_excitatory/run_perturb_array_finetuned_control_to_tbi_overexpress.sh   # induction / overexpress, --array=0-71
# resubmit tasks that hit the wall-time limit:  sbatch --array=<failed ids> <the same script>
# aggregation into per-gene tables (CPU; 8 CPUs, 64 GB, 4 h limit; 13-14 min per arm on the cluster)
\$GF_PYTHON ${S}/03_in_silico_perturbation/ft_excitatory/05_get_perturb_stats_finetuned.py --mode delete --celltype excit_rescue   # rescue / delete (no wrapper)
sbatch ${S}/03_in_silico_perturbation/ft_excitatory/run_05_07_09_finetuned_rescue_overexpress.sh                 # rescue / overexpress: stages 05, 07, 09
sbatch ${S}/03_in_silico_perturbation/ft_excitatory/run_05_get_perturb_stats_finetuned_delete.sh                 # induction / delete
sbatch ${S}/03_in_silico_perturbation/ft_excitatory/run_05_get_perturb_stats_finetuned_overexpress.sh            # induction / overexpress
# The per-gene tables written here are the ones shipped in data/${S}/in_silico_perturbation/.

# ---- then the CPU half, reading the tree just built
bash ${S}/run_section.sh --with-project-tree

# ---- MT model: the six *_RECONSTRUCTED scripts (mt_01 ... mt_06) and figures/mt_cls_umap_RECONSTRUCTED.py
# were rewritten from the Methods and were not used for the published results. They cannot run
# end to end: the MT checkpoint is not available and no script writes the tokenized MT corpora.
# Intended order: mt_01_depth_match_garza.py, mt_02_donor_split_garza.py, (tokenization: missing),
# mt_03_finetune_garza.py, mt_04_build_state_embs.py, mt_05_in_silico_perturbation.py --arm <arm>,
# mt_06_aggregate_isp.py --arm <arm>.
EOF
}

if [ "${PLAN}" = 1 ]; then
  print_cluster_plan
  exit 0
fi

mkdir -p "${LOGS}"
export MPLBACKEND=Agg

SUMMARY=()
NFAIL=0

wanted() {                     # wanted <step name>: true unless --only excludes the step
  [ -z "${ONLY}" ] || [[ "$1" =~ ${ONLY} ]]
}

note() {                       # note <PASS|FAIL|SKIP|SAME|DIFF> <text>
  local line
  line="$(printf '%-4s  %s' "$1" "$2")"
  echo "${line}"
  SUMMARY+=("${line}")
}

# step <log name> <command ...>
# Runs the command, tees stdout+stderr to ${LOGS}/<log name>.log, records PASS/FAIL and the wall time.
step() {
  local name="$1"; shift
  wanted "${name}" || return 0
  echo
  echo "=== ${name}"
  local t0=${SECONDS} rc=0
  "$@" 2>&1 | tee "${LOGS}/${name}.log" || rc=$?
  local dt=$(( SECONDS - t0 ))
  if [ "${rc}" -eq 0 ]; then
    note PASS "${name}  (${dt} s)"
  else
    note FAIL "${name}  (exit ${rc} after ${dt} s, see ${LOGS}/${name}.log)"
    NFAIL=$(( NFAIL + 1 ))
  fi
  return 0
}

# net_step <log name> <command ...>      a step that needs internet access
net_step() {
  local name="$1"; shift
  wanted "${name}" || return 0
  if [ "${OFFLINE}" = 1 ]; then
    note SKIP "${name}  (queries a web service, skipped by --offline)"
    return 0
  fi
  step "${name}" "$@"
}

# skip <log name> <reason>               a step this runner never runs
skip() {
  wanted "$1" || return 0
  note SKIP "$1  ($2)"
}

# tree_step <log name> <reason> <required path> [<required path> ...] -- <command ...>
# A step whose inputs are not shipped in data/. Without --with-project-tree it is skipped with
# the reason. With it, the step runs if every required path exists, else it is skipped naming
# the first one that does not.
tree_step() {
  local name="$1" reason="$2"; shift 2
  local need=()
  while [ $# -gt 0 ] && [ "$1" != "--" ]; do need+=("$1"); shift; done
  shift
  wanted "${name}" || return 0
  if [ "${WITH_TREE}" != 1 ]; then
    note SKIP "${name}  (${reason})"
    return 0
  fi
  local p
  for p in "${need[@]}"; do
    if [ ! -e "${p}" ]; then
      note SKIP "${name}  (${reason}; not found: ${p})"
      return 0
    fi
  done
  step "${name}" "$@"
}

# same <label> <produced file> <shipped file>   byte comparison with the shipped table; never a FAIL
same() {
  wanted "$1" || return 0
  if [ ! -f "$2" ]; then
    note DIFF "$1  (not produced: $2)"
  elif cmp -s "$2" "$3"; then
    note SAME "$1  (byte-identical to data/${3#"${NANOAI_DATA}/"})"
  else
    note DIFF "$1  (differs from data/${3#"${NANOAI_DATA}/"}; expected only when its inputs were regenerated)"
  fi
}

# Locations of inputs that are not shipped (all below the three project roots of gf_env.sh).
DET="${GF_JOINTPG}/results/isp/garza_disease_binary_TBI_to_control/mechano_analysis/layer3/garza_qc4_de_by_celltype.tsv"
MT_TABLE="${GF_JOINTPG}/results/isp/garza_disease_binary_TBI_to_control/stats_delete/by_cell_type/isp_garza_disease_binary_TBI_to_control_delete_excitatory_neuron_unfiltered_all_gene_ranked.tsv"
MT_OLIGO="${GF_JOINTPG}/results/isp/garza_disease_binary_TBI_to_control/stats_delete/by_cell_type/isp_garza_disease_binary_TBI_to_control_delete_oligodendrocyte_filtered_high_confidence_gene_ranked.tsv"
ZS_OLIGO="${GF_WORK}/results/perturb_oligo_paperQC_nGene1000/stats/oligo_paperQC_nGene1000_filtered_high_confidence_gene_ranked.tsv"
H5AD_ALL="${GF_WORK}/h5ad/Garza_GSE209552_all_raw_counts_paperQC_nGene1000_annotated.h5ad"
H5AD_EXC="${GF_WORK}/h5ad/celltype_subsets_paperQC_nGene1000/excitatory_neurons_raw_geneformer_ready.h5ad"
CHUNKS="${GF_WORK}/results/perturb_neurons_finetuned_tbi_to_control/chunks_delete"
EMB_CACHE="${OUT}/umap/excit_ft_cellembs.parquet"
[ -e "${EMB_CACHE}" ] || EMB_CACHE="${GF_WORK}/results/umap/excit_ft_cellembs.parquet"

FT_ISP="${HERE}/03_in_silico_perturbation/ft_excitatory"
P04="${HERE}/04_pathway_enrichment"
P05="${HERE}/05_drug_repurposing"
P06="${HERE}/06_confounders_and_robustness"
FIG="${HERE}/figures"
NO_DET="needs the 20%-detection table garza_qc4_de_by_celltype.tsv, which is not shipped (build_detection_table.py rebuilds it from the stage-01 AnnData)"

echo "repository : ${NANOAI_REPO}"
echo "data       : ${GF_DATA}"
echo "results    : ${OUT}"
echo "work tree  : ${GF_WORK}   (stages 07 and 09 write here)"
echo "python     : $(command -v "${PYTHON}" || echo "${PYTHON} (not found)")"

# ---------------------------------------------------------------- 00-03: the upstream half
skip "upstream_00_to_03" "environment, dataset preparation, fine-tuning, in-silico perturbation and stage 05 need GEO downloads, a GPU and a Slurm cluster; their per-gene tables are shipped in data/. Print the commands with --print-cluster-plan"

# ---------------------------------------------------------------- 03: biological-candidate filter (stage 07)
# CPU-only and reproducible from the shipped per-gene tables: they are first copied into the
# working-tree layout the pipeline stages read (never overwriting a table already there).
step "seed_ft_stats" "${PYTHON}" "${HERE}/gf_paths.py" --seed-ft-stats
# arm : --celltype : --mode
ARMS=("rescue_delete:excit_rescue:delete" "rescue_overexpress:excit_rescue:overexpress"
      "induction_delete:excit:delete" "induction_overexpress:excit:overexpress")
for A in "${ARMS[@]}"; do
  IFS=: read -r ARM CT MODE <<<"${A}"
  step "07_biological_candidate_filter_${ARM}" "${PYTHON}" "${FT_ISP}/07_biological_candidate_filter_finetuned.py" --mode "${MODE}" --celltype "${CT}"
done

# ---------------------------------------------------------------- 04_pathway_enrichment/
for A in "${ARMS[@]}"; do
  IFS=: read -r ARM CT MODE <<<"${A}"
  net_step "09_pathway_enrichment_${ARM}" "${PYTHON}" "${P04}/09_pathway_enrichment_finetuned.py" --mode "${MODE}" --celltype "${CT}"
done
# The detection table of the original run is lost with the MT project; build_detection_table.py
# (not part of the original run) rebuilds its excitatory-neuron rows from the stage-01 AnnData.
if [ -e "${DET}" ]; then
  skip "build_detection_table" "the detection table is already present: ${DET}"
else
  tree_step "build_detection_table" "rebuilds the 20%-detection table; needs the excitatory-neuron AnnData from 01_dataset_preparation, not shipped" \
    "${H5AD_EXC}" -- "${GF_PREP_PYTHON}" "${P04}/build_detection_table.py"
fi
tree_step "excit_pathway_heatmap_enrich" "${NO_DET}; the MT arms also need the MT per-gene tables, which are not available" \
  "${DET}" -- "${PYTHON}" "${P04}/excit_pathway_heatmap_enrich.py"
step "excit_pathway_heatmap_FT" "${PYTHON}" "${P04}/excit_pathway_heatmap_FT.py"
step "excit_pathway_heatmap_MT" "${PYTHON}" "${P04}/excit_pathway_heatmap_MT.py"

# ---------------------------------------------------------------- 05_drug_repurposing/
# union before core, as in the original run: both modes write <model>/all_sig_drug_terms.tsv and
# the shipped copy is the one left by the core run.
for MODE in union core; do
  tree_step "excit_drug_repurpose_FT_${MODE}" "${NO_DET}" \
    "${DET}" -- "${PYTHON}" "${P05}/excit_drug_repurpose.py" FT "${MODE}"
done
for MODE in union core; do
  tree_step "excit_drug_repurpose_MT_${MODE}" "${NO_DET}, and the MT per-gene tables, which are not available" \
    "${DET}" "${MT_TABLE}" -- "${PYTHON}" "${P05}/excit_drug_repurpose.py" MT "${MODE}"
done
step "excit_drug_crossmodel" "${PYTHON}" "${P05}/excit_drug_crossmodel.py"
step "excit_drug_repurpose_figure" "${PYTHON}" "${P05}/excit_drug_repurpose_figure.py"

# ---------------------------------------------------------------- 06_confounders_and_robustness/
tree_step "confound_check_sex_age_drivers" "needs the excitatory-neuron AnnData from 01_dataset_preparation and the per-arm driver lists written by excit_pathway_heatmap_enrich.py, neither shipped" \
  "${H5AD_EXC}" -- "${GF_PREP_PYTHON}" "${P06}/confound_check_sex_age_drivers.py"
net_step "cross_direction_network_analysis" "${PYTHON}" "${P06}/cross_direction_network_analysis.py"
tree_step "compare_singlect_vs_jointpg" "needs the zero-shot single-cell-type tables, the MT per-cell-type tables and the mechano gene list, none shipped and no script in the repository writes the first two" \
  "${ZS_OLIGO}" "${MT_OLIGO}" "${GF_MECHANO}/data/mechano_genes_list.txt" -- "${PYTHON}" "${P06}/compare_singlect_vs_jointpg.py"

# ---------------------------------------------------------------- figures/
skip "mt_cls_umap_RECONSTRUCTED" "MT embedding UMAP: needs a GPU and the MT checkpoint, which is not available; the script is a reconstruction"
tree_step "umap_all_celltypes_expression" "Supp. Fig. 5a: needs the 15,172-nucleus AnnData written by 01_dataset_preparation, not shipped" \
  "${H5AD_ALL}" -- "${GF_PREP_PYTHON}" "${FIG}/umap_all_celltypes_expression.py"
tree_step "umap_excit_ft" "Supp. Fig. 5b: needs a GPU and the fine-tuned FT checkpoint, or the per-nucleus embedding cache excit_ft_cellembs.parquet of the original run; neither shipped" \
  "${EMB_CACHE}" -- "${GF_PYTHON}" "${FIG}/umap_excit_ft.py"
step "make_suppfig5_model_metrics" "${PYTHON}" "${FIG}/make_suppfig5_model_metrics.py"
step "make_fig2_panelB" "${PYTHON}" "${FIG}/make_fig2_panelB.py"
tree_step "heatmap_excit_finetuned_cosine_shifts" "needs the raw perturbation chunk pickles (1.4 GB) and the 20%-detection table, neither shipped" \
  "${CHUNKS}" "${DET}" -- "${GF_PYTHON}" "${FIG}/heatmap_excit_finetuned_cosine_shifts.py"

# ---------------------------------------------------------------- regenerated tables against data/
same "same_pathway_matrix_FT" "${OUT}/pathway_heatmap/excit_pathway_enrichment_matrix_FT.tsv" "${GF_DATA}/pathway_enrichment/excit_pathway_enrichment_matrix_FT.tsv"
same "same_pathway_matrix_MT" "${OUT}/pathway_heatmap/excit_pathway_enrichment_matrix_MT.tsv" "${GF_DATA}/pathway_enrichment/excit_pathway_enrichment_matrix_MT.tsv"
same "same_crossmodel_consensus" "${OUT}/drug_repurpose/crossmodel_consensus.tsv" "${GF_DATA}/drug_signature_enrichment/crossmodel_consensus.tsv"

# ---------------------------------------------------------------- summary
echo
echo "================ ${SECTION}: summary ================"
if [ "${#SUMMARY[@]}" -eq 0 ]; then
  echo "no step matched --only '${ONLY}'"
  exit 2
fi
printf '%s\n' "${SUMMARY[@]}"
echo "results  : ${OUT}"
echo "work tree: ${GF_WORK}/results   (stage 07 and stage 09 outputs)"
echo "logs     : ${LOGS}"
if [ "${NFAIL}" -gt 0 ]; then
  echo "${NFAIL} step(s) FAILED"
  exit 1
fi
echo "all steps that ran PASSED"
