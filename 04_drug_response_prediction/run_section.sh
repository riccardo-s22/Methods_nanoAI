#!/usr/bin/env bash
# Results 4 - drug response prediction (Fig. 4b, 4g; Supp. Figs. 7, 8, 10b-c): rerun everything
# that can be rerun from the shipped data on an ordinary CPU machine.
#
#   bash 04_drug_response_prediction/run_section.sh [--skip-slow] [--sequential] [--quiet] [--with-refit]
#
# Works from any working directory and needs no environment variable.
#   inputs    <repo>/data/04_drug_response_prediction      (override NANOAI_DATA)
#   outputs   <repo>/results/04_drug_response_prediction   (override NANOAI_RESULTS)
#   logs      <outputs>/logs/<step>.log, and the status table <outputs>/logs/run_status.tsv
#
#   --skip-slow    skip the four 20,000-permutation steps (exp3_perm_power.py, its true-log copy and
#                  the SVG re-render of each; about 40 min of CPU each). They run by default.
#   --sequential   run those four slow steps one after the other, in README order. By default, on a
#                  machine with at least 4 CPUs, they are started first and run concurrently in the
#                  background while the other steps proceed (--parallel forces this).
#   --quiet        send script output to the logs only; print just the PASS/FAIL/SKIP lines
#   --with-refit   also attempt the scripts that refit the peaks from the raw spectra, which crash
#                  with NumPy 2.5 (see README, "Known limitations"). With an interpreter that still
#                  has numpy.trapz (PYTHON=/path/to/python) they are run without this flag.
#
# Optional inputs that are not in the repository; the steps needing them print SKIP when absent:
#   EXP3_WB_XLSX   western-blot densitometry workbook corr_wb_swcnt.xlsx (Supp. Fig. 10b-c)
#   EXP3_DIR       a folder holding both the 3 h and the 6 h read (exp3_load.py self-check,
#                  exp3_analysis.py, exp3_build_cache.py)
#
# Exit status is non-zero if any step failed. All steps are attempted regardless.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../shared/env.sh"
source "${HERE}/run_lib.sh"
parse_flags "$@"
status_begin

PY="$PYTHON"
MAKE_SVG="${NANOAI_REPO}/shared/make_svg.py"   # usage: make_svg.py <script> <cwd> <out dir> [ENV=VAL ...]
SVG_OUT="${SECTION_OUT}/svg"

echo "section : ${SECTION_NAME}"
echo "python  : $(command -v "$PY") ($("$PY" -c 'import sys, numpy, scipy, sklearn, pandas, matplotlib; print("Python %s, NumPy %s, SciPy %s, scikit-learn %s, pandas %s, matplotlib %s" % (sys.version.split()[0], numpy.__version__, scipy.__version__, sklearn.__version__, pandas.__version__, matplotlib.__version__))'))"
echo "data    : ${NANOAI_DATA}/${SECTION_NAME}"
echo "results : ${SECTION_OUT}"
echo

REFIT=0
if numpy_has_trapz || [[ "$WITH_REFIT" == 1 ]]; then REFIT=1; fi

slow_main() {
  # Supp. Fig. 7: 20,000 label permutations for each of the two feature sets
  slow_step exp3_perm_power      "$PY" "${HERE}/exp3_perm_power.py"
  # Supp. Fig. 7 as editable SVG panels: the same script re-run under the SVG harness
  slow_step exp3_perm_power__svg "$PY" "$MAKE_SVG" "${HERE}/exp3_perm_power.py" "$SECTION_OUT" "$SVG_OUT" EXP3_TP=3H
}

# In parallel mode the slow steps are started first so that they overlap with everything else.
if [[ "$PARALLEL" == 1 ]]; then slow_main; fi

# ---------------------------------------------------------------- pipeline, in README run order
if have_6h_read; then
  run_step exp3_load "$PY" "${HERE}/exp3_load.py"
else
  skip_step exp3_load "loader self-check: ${NO_6H_REASON}. The loader itself is exercised by every step below"
fi

if [[ "$REFIT" == 1 ]]; then
  run_step exp3_peakfit        "$PY" "${HERE}/exp3_peakfit.py"
  run_step exp3_peakfit_ratios "$PY" "${HERE}/exp3_peakfit_ratios.py"
else
  skip_step exp3_peakfit        "$NO_TRAPZ_REASON"
  skip_step exp3_peakfit_ratios "$NO_TRAPZ_REASON"
fi

if [[ "$REFIT" != 1 ]]; then
  skip_step exp3_build_cache "$NO_TRAPZ_REASON"
elif ! have_6h_read; then
  skip_step exp3_build_cache "$NO_6H_REASON"
else
  run_step exp3_build_cache "$PY" "${HERE}/exp3_build_cache.py"     # writes to results/, never to data/
fi

run_step exp3_export_peakfeats "$PY" "${HERE}/exp3_export_peakfeats.py"

if have_6h_read; then
  run_step exp3_analysis "$PY" "${HERE}/exp3_analysis.py"
else
  skip_step exp3_analysis "$NO_6H_REASON"
fi

# Fig. 4b and Supp. Fig. 8
run_step exp3_drug_axis      "$PY" "${HERE}/exp3_drug_axis.py"
run_step exp3_drug_axis__svg "$PY" "$MAKE_SVG" "${HERE}/exp3_drug_axis.py" "$SECTION_OUT" "$SVG_OUT" EXP3_TP=3H

run_step exp3_drug_effects "$PY" "${HERE}/exp3_drug_effects.py"

# Fig. 4g scores, and the input of figures/plot_wb_swcnt_corr.py
run_step exp3_3H_organoid_score "$PY" "${HERE}/exp3_3H_organoid_score.py"

if [[ "$PARALLEL" != 1 ]]; then slow_main; fi

# ---------------------------------------------------------------- figures/
if have_wb_xlsx; then
  run_step plot_wb_swcnt_corr env EXP3_WB_XLSX="$EXP3_WB_XLSX_PATH" "$PY" "${HERE}/figures/plot_wb_swcnt_corr.py"
else
  skip_step plot_wb_swcnt_corr "Supp. Fig. 10b-c: ${NO_WB_REASON}"
fi

# ---------------------------------------------------------------- supporting_checks/
run_step exp3_3H_pressure "$PY" "${HERE}/supporting_checks/exp3_3H_pressure.py"
run_step exp3_drug_swap   "$PY" "${HERE}/supporting_checks/exp3_drug_swap.py"

# ---------------------------------------------------------------- sensitivity_true_log_ratios/
# Its own runner records its steps in the same status table.
nfail_before=$(grep -c '^FAIL' "$NANOAI_STATUS_FILE" || true)
if ! bash "${HERE}/sensitivity_true_log_ratios/run_all.sh"; then
  nfail_after=$(grep -c '^FAIL' "$NANOAI_STATUS_FILE" || true)
  # a non-zero exit with no failed step recorded means the runner itself broke
  [[ "$nfail_after" -gt "$nfail_before" ]] || _record FAIL "sensitivity_true_log_ratios/run_all.sh" 0 "runner error"
fi

# ---------------------------------------------------------------- compare with the published values
wait   # background slow steps (--parallel) must have written their logs
run_step verify_outputs "$PY" "${HERE}/verify_outputs.py"

status_end || exit 1
