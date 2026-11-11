#!/usr/bin/env bash
# Re-run of the chirality-ratio analyses of the drug experiment with true log10 ratios (flo = 0).
# The published (floored) scripts and their outputs are untouched.
#
#   bash 04_drug_response_prediction/sensitivity_true_log_ratios/run_all.sh [--skip-slow] [--sequential] [--quiet]
#
# Works from any working directory and needs no environment variable. Called by ../run_section.sh,
# and can be run on its own.
#   inputs    <repo>/data/04_drug_response_prediction                                (override NANOAI_DATA)
#   outputs   <repo>/results/04_drug_response_prediction/sensitivity_true_log_ratios (override NANOAI_RESULTS)
#   logs      <repo>/results/04_drug_response_prediction/logs/<step>.log
#
# The blast-classification half of the original re-run (exp2_chirality_classify_truelog.py) lives in
# 01_swcnt_blast_classification/sensitivity_true_log_ratios/ and is run by that section's runner.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../../shared/env.sh"
source "${HERE}/../run_lib.sh"
parse_flags "$@"
status_begin

PY="$PYTHON"
MAKE_SVG="${NANOAI_REPO}/shared/make_svg.py"
TL_OUT="${SECTION_OUT}/sensitivity_true_log_ratios"

slow_truelog() {
  # 20,000 label permutations per feature set, then the same script under the SVG harness
  slow_step exp3_perm_power_truelog      "$PY" "${HERE}/exp3_perm_power_truelog.py"
  slow_step exp3_perm_power_truelog__svg "$PY" "$MAKE_SVG" "${HERE}/exp3_perm_power_truelog.py" "$SECTION_OUT" "${TL_OUT}/svg" EXP3_TP=3H
}

if [[ "$PARALLEL" == 1 ]]; then slow_truelog; fi

run_step exp3_3H_organoid_score_truelog "$PY" "${HERE}/exp3_3H_organoid_score_truelog.py"

run_step exp3_drug_axis_truelog      "$PY" "${HERE}/exp3_drug_axis_truelog.py"
run_step exp3_drug_axis_truelog__svg "$PY" "$MAKE_SVG" "${HERE}/exp3_drug_axis_truelog.py" "$SECTION_OUT" "${TL_OUT}/svg" EXP3_TP=3H

# true log ratios (FLO defaults to 0.0), then the control run with the published floor, which must
# give back the published 0.613 / 0.938 / 0.920
run_step exp3_peakfit_ratios_truelog            "$PY" "${HERE}/exp3_peakfit_ratios_truelog.py"
run_step exp3_peakfit_ratios_truelog__FLO1e-6   env FLO=1e-6 EXP3_OUT="${TL_OUT}/checkfloor_FLO1e-6" "$PY" "${HERE}/exp3_peakfit_ratios_truelog.py"

if have_wb_xlsx; then
  run_step plot_wb_swcnt_corr_truelog env EXP3_WB_XLSX="$EXP3_WB_XLSX_PATH" "$PY" "${HERE}/plot_wb_swcnt_corr_truelog.py"
else
  skip_step plot_wb_swcnt_corr_truelog "$NO_WB_REASON"
fi

if [[ "$PARALLEL" != 1 ]]; then slow_truelog; fi

status_end || exit 1
