#!/usr/bin/env bash
# Results 1 - SWCNT blast classification (Fig. 1, Supp. Figs. 1-4).
# Re-runs every script of this section from the data shipped in data/, in the order of README.md.
#
#     bash 01_swcnt_blast_classification/run_section.sh [--skip-slow] [--skip-svg] [--only REGEX]
#
# Works from any working directory. Nothing is read from or written to the current directory:
#   inputs    <repo>/data/01_swcnt_blast_classification      (NANOAI_DATA overrides <repo>/data)
#   outputs   <repo>/results/01_swcnt_blast_classification   (NANOAI_RESULTS overrides <repo>/results)
#   logs      <repo>/results/01_swcnt_blast_classification/logs/<step>.log
#
#   --skip-slow   skip the steps that take minutes each (permutation tests, classifier comparison,
#                 Supp. Fig. 2 computation, true-log-ratio null). Also: SKIP_SLOW=1.
#   --skip-svg    skip the SVG re-rendering of the published panels through shared/make_svg.py
#                 (it runs each figure script a second time). Also: SKIP_SVG=1.
#   --only REGEX  run only the steps whose name matches the extended regular expression, e.g.
#                 --only 'exp2_mixes'  or  --only '_24h$'  or  --only 'chirality|gaussian'.
#                 Step names are the names in the PASS/FAIL lines and of the log files.
#
# Every step prints one PASS / FAIL / SKIP line; all of them are repeated in a summary at the end.
# The exit status is non-zero if any step failed. A failing step does not stop the run.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../shared/env.sh
source "${HERE}/../shared/env.sh"

SECTION="01_swcnt_blast_classification"
DATA="${NANOAI_DATA}/${SECTION}"
OUT="${NANOAI_RESULTS}/${SECTION}"
LOGS="${OUT}/logs"
SVG="${NANOAI_REPO}/shared/make_svg.py"
SKIP_SLOW="${SKIP_SLOW:-0}"
SKIP_SVG="${SKIP_SVG:-0}"
ONLY=""

while [ $# -gt 0 ]; do
  case "$1" in
    --skip-slow) SKIP_SLOW=1 ;;
    --skip-svg)  SKIP_SVG=1 ;;
    --only)      [ $# -ge 2 ] || { echo "--only needs a regular expression" >&2; exit 2; }
                 ONLY="$2"; shift ;;
    --only=*)    ONLY="${1#--only=}" ;;
    -h|--help)   sed -n '2,21p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1 (use --help)" >&2; exit 2 ;;
  esac
  shift
done

mkdir -p "${LOGS}"
export MPLBACKEND=Agg
# One BLAS/OpenMP thread per process. Every script but exp2_analysis.py pins this itself; on the
# small matrices of this section a multi-threaded BLAS is many times slower, not faster, and the
# printed statistics are identical either way. Set the variables yourself to override.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export MIX_INVERTED=1          # column 10 = 25% blast, column 12 = 75% blast (the published labelling)
unset EXP2_OUT                 # outputs of this run always go under ${OUT}

SUMMARY=()
NFAIL=0

wanted() {                     # wanted <step name>: true unless --only excludes the step
  [ -z "${ONLY}" ] || [[ "$1" =~ ${ONLY} ]]
}

note() {                       # note <PASS|FAIL|SKIP> <text>
  local line
  line="$(printf '%-4s  %s' "$1" "$2")"
  echo "${line}"
  SUMMARY+=("${line}")
}

# step <log name> <fast|slow> <command ...>
# Runs the command, tees stdout+stderr to ${LOGS}/<log name>.log, records PASS/FAIL and the wall time.
step() {
  local name="$1" speed="$2"; shift 2
  wanted "${name}" || return 0
  if [ "${speed}" = slow ] && [ "${SKIP_SLOW}" = 1 ]; then
    note SKIP "${name}  (slow step, skipped by --skip-slow)"
    return 0
  fi
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

# svg <log name> <fast|slow> <script> <out subdir> [ENV=VAL ...]   SVG re-render of a figure script
svg() {
  local name="$1" speed="$2" script="$3" sub="$4"; shift 4
  wanted "${name}" || return 0
  if [ "${SKIP_SVG}" = 1 ]; then
    note SKIP "${name}  (SVG re-render, skipped by --skip-svg)"
    return 0
  fi
  step "${name}" "${speed}" "${PYTHON}" "${SVG}" "${script}" "${SECTION}/svg/${sub}" "$@"
}

spectra() { echo "${DATA}/spectra/$1"; }

echo "repository : ${NANOAI_REPO}"
echo "data       : ${DATA}"
echo "results    : ${OUT}"
echo "python     : $(command -v "${PYTHON}" || echo "${PYTHON} (not found)")"

# ---------------------------------------------------------------- top level: the analysis pipeline
# One run per read of the plate. EXP2_DIR selects the spectra, EXP2_TAG the label in titles and file names.
for TP in 2h 24h; do
  export EXP2_TAG="${TP}" EXP2_DIR; EXP2_DIR="$(spectra "${TP}")"
  step "exp2_load_${TP}"               fast "${PYTHON}" "${HERE}/exp2_load.py"
  step "exp2_analysis_${TP}"           slow "${PYTHON}" "${HERE}/exp2_analysis.py"
  step "exp2_analysis2_${TP}"          slow "${PYTHON}" "${HERE}/exp2_analysis2.py"
  step "exp2_classifier_compare_${TP}" slow "${PYTHON}" "${HERE}/exp2_classifier_compare.py"
  step "exp2_dose_pca_figs_${TP}"      fast "${PYTHON}" "${HERE}/exp2_dose_pca_figs.py"
  step "exp2_mixes_${TP}"              fast "${PYTHON}" "${HERE}/exp2_mixes.py"
  step "exp2_mix_intermediate_${TP}"   fast "${PYTHON}" "${HERE}/exp2_mix_intermediate.py"
  step "exp2_mix_features_${TP}"       fast "${PYTHON}" "${HERE}/exp2_mix_features.py"
done
unset EXP2_TAG EXP2_DIR
step "exp2_2h_vs_24h" fast "${PYTHON}" "${HERE}/exp2_2h_vs_24h.py"      # reads both reads itself

# ---------------------------------------------------------------- SVG panels (Fig. 1, Supp. Fig. 1)
for TP in 2h 24h; do
  E=("EXP2_TAG=${TP}" "EXP2_DIR=$(spectra "${TP}")" "MIX_INVERTED=1")
  svg "svg_exp2_analysis2_${TP}"          slow "${HERE}/exp2_analysis2.py"          "exp2_${TP}" "${E[@]}"
  svg "svg_exp2_classifier_compare_${TP}" slow "${HERE}/exp2_classifier_compare.py" "exp2_${TP}" "${E[@]}"
  svg "svg_exp2_dose_pca_figs_${TP}"      fast "${HERE}/exp2_dose_pca_figs.py"      "exp2_${TP}" "${E[@]}"
  svg "svg_exp2_mixes_${TP}"              fast "${HERE}/exp2_mixes.py"              "exp2_${TP}" "${E[@]}"
  svg "svg_exp2_mix_intermediate_${TP}"   fast "${HERE}/exp2_mix_intermediate.py"   "exp2_${TP}" "${E[@]}"
done

# ---------------------------------------------------------------- figures/
# Supp. Fig. 2: compute (slow) then draw. The figure step needs the JSON the compute step writes.
step "make_supp_auc_compute" slow "${PYTHON}" "${HERE}/figures/make_supp_auc_compute.py"
if [ -f "${OUT}/supp_fig2_auc/supp_auc_results.json" ]; then
  step "make_supp_auc_figure" fast "${PYTHON}" "${HERE}/figures/make_supp_auc_figure.py"
elif wanted "make_supp_auc_figure"; then
  note SKIP "make_supp_auc_figure  (needs ${OUT}/supp_fig2_auc/supp_auc_results.json from make_supp_auc_compute)"
fi

# Supp. Fig. 3 and the Gaussian version of Supp. Fig. 4 need the cached split-Gaussian fits of the
# sham and blast wells, exp2_{2h,24h}_peakfeats.npz. They are looked for in data/ (or in the folder
# named by EXP2_PEAKFEATS_DIR); if they are not there the steps are skipped, not failed.
PEAK_DIR="${EXP2_PEAKFEATS_DIR:-${DATA}/chirality_features}"
HAVE_PEAK=1
for TP in 2h 24h; do
  [ -f "${PEAK_DIR}/exp2_${TP}_peakfeats.npz" ] || HAVE_PEAK=0
done
NEED_PEAK="input not in the repository: exp2_2h_peakfeats.npz and exp2_24h_peakfeats.npz (sham/blast split-Gaussian fits) not found in ${PEAK_DIR}; set EXP2_PEAKFEATS_DIR to the folder holding them"

# peak_step <log name> <fast|slow> <command ...>   a step that needs the two peakfeats files
peak_step() {
  local name="$1"
  if [ "${HAVE_PEAK}" = 1 ]; then
    step "$@"
  elif wanted "${name}"; then
    note SKIP "${name}  (${NEED_PEAK})"
  fi
}

peak_step "exp2_chirality_classify" fast "${PYTHON}" "${HERE}/figures/exp2_chirality_classify.py"
step "make_suppfig4_dilution_readouts" fast "${PYTHON}" "${HERE}/figures/make_suppfig4_dilution_readouts.py"
peak_step "make_suppfig4_gaussian_readouts" fast "${PYTHON}" "${HERE}/figures/make_suppfig4_gaussian_readouts.py"
step "exp2_pca_quant" fast "${PYTHON}" "${HERE}/figures/exp2_pca_quant.py"

# ---------------------------------------------------------------- supporting_checks/
step "check_simple_descriptors" fast "${PYTHON}" "${HERE}/supporting_checks/check_simple_descriptors.py"
step "check_permutation_power"  fast "${PYTHON}" "${HERE}/supporting_checks/check_permutation_power.py"

# ---------------------------------------------------------------- sensitivity_true_log_ratios/
# Slow on the first run only: the 5,000-permutation null of the true log ratios is not shipped,
# so it is computed (8 workers) and cached under ${OUT}/sensitivity_true_log_ratios/.
peak_step "exp2_chirality_classify_truelog" slow "${PYTHON}" "${HERE}/sensitivity_true_log_ratios/exp2_chirality_classify_truelog.py"
peak_step "make_suppfig4_gaussian_readouts_truelog" fast env FEAT=truelog "${PYTHON}" "${HERE}/figures/make_suppfig4_gaussian_readouts.py"

# ---------------------------------------------------------------- summary
echo
echo "================ ${SECTION}: summary ================"
if [ "${#SUMMARY[@]}" -eq 0 ]; then
  echo "no step matched --only '${ONLY}'"
  exit 2
fi
printf '%s\n' "${SUMMARY[@]}"
echo "results: ${OUT}"
echo "logs   : ${LOGS}"
if [ "${NFAIL}" -gt 0 ]; then
  echo "${NFAIL} step(s) FAILED"
  exit 1
fi
echo "all steps that ran PASSED"
