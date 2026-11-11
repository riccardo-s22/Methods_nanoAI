#!/usr/bin/env bash
# Rerun every analysis of the paper that can be rerun from this repository.
#
#     bash reproduce.sh                  all four sections, in the order 01 02 03 04
#     bash reproduce.sh 01 04            a subset; sections are given as 01..04 (or 1..4, or the
#                                        directory name) and always run in ascending order
#     bash reproduce.sh --list [01 ...]  print every step of every section with its class and exit;
#                                        nothing is run
#
# It calls <section>/run_section.sh for each section, then compares the regenerated tables with the
# reference tables shipped in data/ (repro/check_results.py), and ends with one table:
# section, steps passed / failed / skipped, wall time.
#
# Options:
#   --skip-slow   forwarded to sections 01 and 04: skip the permutation-heavy steps
#                 (full run about 1 h 50 min, with --skip-slow about 15 min)
#   --skip-svg    forwarded to section 01: skip the SVG re-rendering of the figure panels
#   --offline     no internet: forwarded to section 02 (skips its Enrichr and STRING steps);
#                 section 03 is not run, because it has no offline mode
#   --no-check    do not run repro/check_results.py at the end
#   -h, --help    this text
#
# Works from any working directory. Nothing is read from or written to the current directory:
#   inputs    <repo>/data       read-only        (override with NANOAI_DATA)
#   outputs   <repo>/results    git-ignored      (override with NANOAI_RESULTS)
#             section 02 also writes its stage 07 / 09 tables to <repo>/work (override with GF_WORK)
#   logs      results/<section>/logs/<step>.log  one per step
#             results/reproduce/<section>.log    the full output of each section runner
#             results/reproduce/summary.tsv      the table printed at the end
#             results/reproduce/check_results.txt
#
# Interpreter: ${PYTHON}, default python3 (see envs/). Environment variables understood by the
# section runners (EXP2_PEAKFEATS_DIR, EXP3_WB_XLSX, EXP3_DIR, GF_WORK, ...) are passed through.
#
# Exit status: 0 if every step that ran passed and no regenerated table differs from data/;
# 1 if any step failed or a table differs; 2 for a usage or setup error.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=shared/env.sh
source "${ROOT}/shared/env.sh"

ALL=(01_swcnt_blast_classification 02_geneformer_in_silico_screen
     03_protein_corona_proteomics 04_drug_response_prediction)
STEPS="${ROOT}/repro/steps.tsv"
TOP="${NANOAI_RESULTS}/reproduce"

usage() { awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "${BASH_SOURCE[0]}"; }

resolve() {                    # resolve <01|1|directory name>: prints the section directory
  local s
  for s in "${ALL[@]}"; do
    if [ "$1" = "$s" ] || [ "$1" = "${s%%_*}" ] || [ "0$1" = "${s%%_*}" ]; then
      echo "$s"; return 0
    fi
  done
  return 1
}

fmt_time() {                   # fmt_time <seconds>
  local t="$1"
  if [ "$t" -ge 3600 ]; then printf '%dh%02dm%02ds' $((t / 3600)) $((t % 3600 / 60)) $((t % 60))
  elif [ "$t" -ge 60 ]; then printf '%dm%02ds' $((t / 60)) $((t % 60))
  else printf '%ds' "$t"; fi
}

list_steps() {                 # list_steps <section> ...   from repro/steps.tsv, runs nothing
  local s
  echo "Steps of each section runner, in run order. Nothing is run by --list."
  echo "class:  runs locally             runs from a fresh clone on a CPU machine"
  echo "        needs unshipped input    prints SKIP unless the named input is supplied"
  echo "        needs refit environment  prints SKIP unless PYTHON has numpy.trapz (envs/refit.yml)"
  echo "        needs cluster            never run by the runners (GPU, Slurm, large downloads)"
  echo "[slow]  skipped by --skip-slow"
  for s in "$@"; do
    echo
    echo "== ${s}    (bash ${s}/run_section.sh)"
    awk -F'\t' -v s="$s" '
      NR > 1 && $1 == s {
        printf "  %-24s %-6s %-52s %s\n", $3, ($4 == "slow" ? "[slow]" : ""), $2, $5; n[$3]++
      }
      END {
        printf "  --"
        for (c in n) printf "  %s: %d", c, n[c]
        printf "\n"
      }' "${STEPS}"
  done
}

count() {                      # count <PASS|FAIL|SKIP|SAME|DIFF> <runner log>: lines of the final summary
  awk -v k="$1  " '/^=+ .*: summary =+$/ { on = 1; next } on && index($0, k) == 1 { n++ } END { print n + 0 }' "$2"
}

check_manifest() {             # warn when the steps a runner reported are not those of repro/steps.tsv
  local dir="$1" log="$2" diff
  diff="$(comm -3 \
    <(awk '/^=+ .*: summary =+$/ { on = 1; next } on && /^(PASS|FAIL|SKIP)  / { print $2 }' "$log" | sort -u) \
    <(awk -F'\t' -v s="$dir" 'NR > 1 && $1 == s && $2 !~ /^upstream:/ { print $2 }' "${STEPS}" | sort -u))"
  if [ -n "${diff}" ]; then
    echo "note: the steps reported by ${dir}/run_section.sh differ from repro/steps.tsv (used by --list):"
    echo "${diff}" | sed 's/^/        /'
  fi
}

main() {
  local LIST=0 CHECK=1 SKIP_SLOW=0 SKIP_SVG=0 OFFLINE=0
  local want=() sections=() a s dir

  while [ $# -gt 0 ]; do
    case "$1" in
      --list)      LIST=1 ;;
      --skip-slow) SKIP_SLOW=1 ;;
      --skip-svg)  SKIP_SVG=1 ;;
      --offline)   OFFLINE=1 ;;
      --no-check)  CHECK=0 ;;
      -h|--help)   usage; return 0 ;;
      -*)          echo "unknown option: $1 (use --help)" >&2; return 2 ;;
      *)           dir="$(resolve "$1")" || { echo "unknown section: $1 (use 01, 02, 03 or 04)" >&2; return 2; }
                   want+=("${dir}") ;;
    esac
    shift
  done
  if [ "${#want[@]}" -eq 0 ]; then
    sections=("${ALL[@]}")
  else
    for s in "${ALL[@]}"; do               # ascending order, each section once
      for a in "${want[@]}"; do
        if [ "$a" = "$s" ]; then sections+=("$s"); break; fi
      done
    done
  fi

  if [ "${LIST}" = 1 ]; then
    list_steps "${sections[@]}"
    return 0
  fi

  if ! "${PYTHON}" -c 'import numpy, scipy, sklearn, pandas, matplotlib' 2>/dev/null; then
    echo "The interpreter '${PYTHON}' cannot import numpy, scipy, sklearn, pandas and matplotlib." >&2
    echo "Create the environment in envs/ (see README.md, Reproducing the paper), activate it," >&2
    echo "or set PYTHON=/path/to/python." >&2
    return 2
  fi

  mkdir -p "${TOP}"
  echo "repository : ${NANOAI_REPO}"
  echo "data       : ${NANOAI_DATA}"
  echo "results    : ${NANOAI_RESULTS}"
  echo "python     : $(command -v "${PYTHON}")  ($("${PYTHON}" -c 'import sys, numpy, scipy, sklearn, pandas, matplotlib; print("Python %s, NumPy %s, SciPy %s, scikit-learn %s, pandas %s, matplotlib %s" % (sys.version.split()[0], numpy.__version__, scipy.__version__, sklearn.__version__, pandas.__version__, matplotlib.__version__))'))"
  echo "sections   : ${sections[*]}"

  local rows=() bad=0 t_all=${SECONDS}
  local args log t0 dt rc np nf ns extra status
  for dir in "${sections[@]}"; do
    args=()
    case "${dir}" in
      01_*) [ "${SKIP_SLOW}" = 1 ] && args+=(--skip-slow); [ "${SKIP_SVG}" = 1 ] && args+=(--skip-svg) ;;
      02_*) [ "${OFFLINE}" = 1 ] && args+=(--offline) ;;
      04_*) [ "${SKIP_SLOW}" = 1 ] && args+=(--skip-slow) ;;
    esac
    echo
    echo "################################################################ ${dir}"
    if [ "${dir}" = 03_protein_corona_proteomics ] && [ "${OFFLINE}" = 1 ]; then
      echo "SKIP  ${dir}: not run with --offline (ms_analysis.py needs UniProt, g:Profiler and Enrichr)"
      rows+=("${dir}"$'\t'"0"$'\t'"0"$'\t'"2"$'\t'"0s"$'\t'"not run (--offline)")
      continue
    fi
    log="${TOP}/${dir}.log"
    t0=${SECONDS}
    bash "${ROOT}/${dir}/run_section.sh" ${args[@]+"${args[@]}"} 2>&1 | tee "${log}"
    rc=${PIPESTATUS[0]}
    dt=$(( SECONDS - t0 ))
    extra=""
    if grep -Eq '^=+ .*: summary =+$' "${log}"; then
      np="$(count PASS "${log}")"; nf="$(count FAIL "${log}")"; ns="$(count SKIP "${log}")"
      if [ "${dir}" = 02_geneformer_in_silico_screen ]; then
        extra="tables identical to data/: $(count SAME "${log}"), different: $(count DIFF "${log}")"
      fi
      check_manifest "${dir}" "${log}"
    else
      np="?"; nf="?"; ns="?"
      extra="runner stopped before its summary"
    fi
    status="ok"
    if [ "${rc}" -ne 0 ] || [ "${nf}" != 0 ]; then
      status="FAILED (exit ${rc})"; bad=1
    fi
    [ -n "${extra}" ] && status="${status}; ${extra}"
    rows+=("${dir}"$'\t'"${np}"$'\t'"${nf}"$'\t'"${ns}"$'\t'"$(fmt_time "${dt}")"$'\t'"${status}")
  done

  local check_line="not run (--no-check)"
  if [ "${CHECK}" = 1 ]; then
    echo
    echo "################################################################ repro/check_results.py"
    "${PYTHON}" "${ROOT}/repro/check_results.py" 2>&1 | tee "${TOP}/check_results.txt"
    rc=${PIPESTATUS[0]}
    check_line="$(tail -n 1 "${TOP}/check_results.txt")"
    if [ "${rc}" -ne 0 ]; then bad=1; fi
  fi

  {
    printf 'section\tpassed\tfailed\tskipped\twall_time\tstatus\n'
    printf '%s\n' "${rows[@]}"
  } > "${TOP}/summary.tsv"

  echo
  echo "================================ reproduce.sh: summary ================================"
  awk -F'\t' '{ printf "%-32s %7s %7s %8s %10s   %s\n", $1, $2, $3, $4, $5, $6 }' "${TOP}/summary.tsv"
  echo "total wall time : $(fmt_time $(( SECONDS - t_all )))"
  echo "check_results   : ${check_line}"
  echo "outputs         : ${NANOAI_RESULTS}/<section>/    logs: ${NANOAI_RESULTS}/<section>/logs/"
  echo "this summary    : ${TOP}/summary.tsv"
  echo "expected values : ${ROOT}/docs/EXPECTED_RESULTS.md"
  echo "SKIP lines are expected on a fresh clone: see 'bash reproduce.sh --list' for the reason of each."
  if [ "${bad}" -ne 0 ]; then
    echo "RESULT: at least one step FAILED or one regenerated table differs from data/"
    return 1
  fi
  echo "RESULT: every step that ran PASSED"
  return 0
}

main "$@"; exit $?
