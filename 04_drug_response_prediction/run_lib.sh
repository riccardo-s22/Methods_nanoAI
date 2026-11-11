#!/usr/bin/env bash
# Step helpers shared by run_section.sh and sensitivity_true_log_ratios/run_all.sh.
# Source it after shared/env.sh; do not execute it.
#
#   run_step  <name> <command...>   run, tee stdout+stderr to logs/<name>.log, print PASS/FAIL
#   slow_step <name> <command...>   same, but skipped by --skip-slow, and run in the background
#                                   (concurrently with the other slow steps) unless --sequential
#   skip_step <name> <reason>       print one SKIP line
#   status_begin / status_end       open and close the status table of one invocation
#
# Every step appends one line to logs/run_status.tsv: STATUS <tab> name <tab> seconds <tab> note.

SECTION_NAME="04_drug_response_prediction"
SECTION_DIR="${NANOAI_REPO}/${SECTION_NAME}"
SECTION_OUT="${NANOAI_RESULTS}/${SECTION_NAME}"
LOGDIR="${SECTION_OUT}/logs"
mkdir -p "$LOGDIR"

export MPLBACKEND=Agg                 # no display needed
export PYTHONDONTWRITEBYTECODE=1      # keep the source tree free of __pycache__
export EXP3_TP="${EXP3_TP:-3H}"       # the read analysed in the paper
# One BLAS/OpenMP thread per process. The analysis scripts pin this themselves at their top, but
# under shared/make_svg.py NumPy is already imported by then, so the SVG re-renders would run a
# multi-threaded BLAS: much slower on these small matrices, and it oversubscribes the CPUs when
# the four slow steps run concurrently. Set the variables yourself to override.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

SKIP_SLOW="${SKIP_SLOW:-0}"           # 1 = do not run the 20,000-permutation steps
PARALLEL="${PARALLEL:-auto}"          # 1 = run the slow steps concurrently in the background, 0 = one after
                                      # the other; auto = 1 on a machine with at least 4 CPUs
QUIET="${QUIET:-0}"                   # 1 = script output goes to the log only
WITH_REFIT="${WITH_REFIT:-0}"         # 1 = attempt the refit-from-raw scripts even where they are known to crash

parse_flags() {
  local a
  for a in "$@"; do
    case "$a" in
      --skip-slow)  SKIP_SLOW=1 ;;
      --parallel)   PARALLEL=1 ;;
      --sequential) PARALLEL=0 ;;
      --quiet)      QUIET=1 ;;
      --with-refit) WITH_REFIT=1 ;;
      -h|--help)    sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
      *) echo "unknown option: $a (try --help)" >&2; exit 2 ;;
    esac
  done
  if [[ "$PARALLEL" == auto ]]; then
    # nproc reports OMP_NUM_THREADS when that is set (it is, above), so ask without it
    if [[ "$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc 2>/dev/null || echo 1)" -ge 4 ]]; then PARALLEL=1; else PARALLEL=0; fi
  fi
  export SKIP_SLOW PARALLEL QUIET WITH_REFIT
}

status_begin() {
  # The outermost runner owns the status table; a runner called from it appends to the same one.
  if [[ -z "${NANOAI_STATUS_FILE:-}" ]]; then
    export NANOAI_STATUS_FILE="${LOGDIR}/run_status.tsv"
    : > "$NANOAI_STATUS_FILE"
    STATUS_OWNER=1
  else
    STATUS_OWNER=0
  fi
  STATUS_START_LINE=$(wc -l < "$NANOAI_STATUS_FILE")
}

_record() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" >> "$NANOAI_STATUS_FILE"; }

_run() {  # _run <name> <to_terminal 0|1> <command...>
  local name="$1" show="$2"; shift 2
  local log="${LOGDIR}/${name}.log" t0=$SECONDS rc=0
  if [[ "$show" == 1 ]]; then
    "$@" 2>&1 | tee "$log" || rc=$?
  else
    "$@" > "$log" 2>&1 || rc=$?
  fi
  local dt=$((SECONDS - t0))
  if [[ $rc -eq 0 ]]; then
    echo "PASS  ${name}  (${dt}s)  log: ${log}"; _record PASS "$name" "$dt" ""
  else
    echo "FAIL  ${name}  (${dt}s, exit ${rc})  log: ${log}"; _record FAIL "$name" "$dt" "exit ${rc}"
  fi
  return 0
}

run_step() {
  local name="$1"; shift
  echo "----- ${name}"
  if [[ "$QUIET" == 1 ]]; then _run "$name" 0 "$@"; else _run "$name" 1 "$@"; fi
}

skip_step() {
  echo "SKIP  $1  ($2)"; _record SKIP "$1" 0 "$2"
}

slow_step() {
  local name="$1"; shift
  if [[ "$SKIP_SLOW" == 1 ]]; then
    skip_step "$name" "slow step, skipped by --skip-slow"
  elif [[ "$PARALLEL" == 1 ]]; then
    echo "START ${name}  (background, output in ${LOGDIR}/${name}.log)"
    _run "$name" 0 "$@" &
  else
    run_step "$name" "$@"
  fi
}

status_end() {
  wait   # background slow steps, if any
  local mine nfail
  mine=$(tail -n +"$((STATUS_START_LINE + 1))" "$NANOAI_STATUS_FILE")
  nfail=$(printf '%s\n' "$mine" | grep -c '^FAIL' || true)
  if [[ "$STATUS_OWNER" == 1 ]]; then
    echo
    echo "================ ${SECTION_NAME}: summary ================"
    awk -F'\t' '{ printf "%-5s %-46s %6ss  %s\n", $1, $2, $3, $4 }' "$NANOAI_STATUS_FILE"
    echo "PASS $(grep -c '^PASS' "$NANOAI_STATUS_FILE" || true)   FAIL $(grep -c '^FAIL' "$NANOAI_STATUS_FILE" || true)   SKIP $(grep -c '^SKIP' "$NANOAI_STATUS_FILE" || true)"
    echo "outputs: ${SECTION_OUT}"
    echo "logs:    ${LOGDIR}   (status table: ${NANOAI_STATUS_FILE})"
  fi
  [[ "$nfail" -eq 0 ]]
}

# --- availability of optional inputs -------------------------------------------------------------

EXP3_SPECTRA_DIR="${EXP3_DIR:-${NANOAI_DATA}/${SECTION_NAME}/spectra}"
EXP3_WB_XLSX_PATH="${EXP3_WB_XLSX:-${NANOAI_DATA}/${SECTION_NAME}/western_blot/corr_wb_swcnt.xlsx}"

have_6h_read() {  # the 6 h read of the same plates: 5 SHAM + 5 TBI files, not shipped in data/
  local n
  n=$(find "$EXP3_SPECTRA_DIR" -maxdepth 1 -name '*6H__*.txt' 2>/dev/null | wc -l)
  [[ "$n" -ge 10 ]]
}

have_wb_xlsx() { [[ -f "$EXP3_WB_XLSX_PATH" ]]; }

numpy_has_trapz() { "$PYTHON" -c 'import numpy, sys; sys.exit(0 if hasattr(numpy, "trapz") else 1)'; }

NO_6H_REASON="needs the 6 h read (SHAM_EX3_6H__*.txt, TBI_EX3_6H__*.txt), which is not in the repository; point EXP3_DIR at a folder holding both reads to run it"
NO_WB_REASON="needs the western-blot densitometry workbook corr_wb_swcnt.xlsx, which is not in the repository; set EXP3_WB_XLSX to run it"
NO_TRAPZ_REASON="refit from raw spectra crashes with the installed NumPy (no numpy.trapz, used by the fit-failure fallback); the pipeline reads the shipped cache instead. It runs with SciPy 1.15 and NumPy 2.2 (set PYTHON to such an interpreter); --with-refit attempts it anyway"
