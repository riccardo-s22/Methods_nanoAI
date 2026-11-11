#!/usr/bin/env bash
# Results 3 - SWCNT protein corona proteomics (Fig. 3).
# Re-runs ms_analysis.py from the Spectronaut report shipped in data/.
#
#     bash 03_protein_corona_proteomics/run_section.sh
#
# Works from any working directory. Nothing is read from or written to the current directory:
#   input     <repo>/data/03_protein_corona_proteomics/Riccado_Report.csv   (NANOAI_DATA overrides <repo>/data)
#   outputs   <repo>/results/03_protein_corona_proteomics                    (NANOAI_RESULTS overrides <repo>/results)
#   logs      <repo>/results/03_protein_corona_proteomics/logs/<step>.log
#
# ms_analysis.py NEEDS INTERNET ACCESS. It queries three public services while it runs:
#   rest.uniprot.org        UniProt REST API         accession -> primary gene symbol
#   biit.cs.ut.ee           g:Profiler REST API      g:Convert (Ensembl IDs) and g:GOSt (enrichment)
#   maayanlab.cloud         Enrichr                  gene-set libraries for the pre-ranked GSEA
#                                                    (downloaded by gseapy, then cached in ~/.cache/gseapy)
# No gene-set or annotation file is shipped: these services are the only source.
#
# Every step prints one PASS / FAIL line; all of them are repeated in a summary at the end.
# The exit status is non-zero if any step failed. A failing step does not stop the run.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../shared/env.sh
source "${HERE}/../shared/env.sh"

SECTION="03_protein_corona_proteomics"
DATA="${NANOAI_DATA}/${SECTION}"
OUT="${NANOAI_RESULTS}/${SECTION}"
LOGS="${OUT}/logs"

for arg in "$@"; do
  case "${arg}" in
    -h|--help) sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: ${arg} (use --help)" >&2; exit 2 ;;
  esac
done

mkdir -p "${LOGS}"
export MPLBACKEND=Agg

SUMMARY=()
NFAIL=0

note() {                       # note <PASS|FAIL|SKIP> <text>
  local line
  line="$(printf '%-4s  %s' "$1" "$2")"
  echo "${line}"
  SUMMARY+=("${line}")
}

# step <log name> <command ...>
# Runs the command, tees stdout+stderr to ${LOGS}/<log name>.log, records PASS/FAIL and the wall time.
step() {
  local name="$1"; shift
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

# Checks what ms_analysis.py needs before it starts: the input report, the Python packages it
# imports, and that the three services answer. Fails with the reason instead of a traceback.
preflight() {
  "${PYTHON}" - "${DATA}/Riccado_Report.csv" <<'PY'
import importlib, sys
ok = True
import os
if os.path.isfile(sys.argv[1]):
    print("input    ok       ", sys.argv[1])
else:
    print("input    MISSING  ", sys.argv[1]); ok = False
for mod in ("numpy", "pandas", "scipy", "statsmodels", "matplotlib", "requests",
            "adjustText", "gseapy", "openpyxl"):
    try:
        m = importlib.import_module(mod)
        print("package  ok        %-12s %s" % (mod, getattr(m, "__version__", "")))
    except Exception as ex:
        print("package  MISSING   %-12s (%s)" % (mod, ex)); ok = False
try:
    import requests
    for name, url in (("UniProt REST", "https://rest.uniprot.org/"),
                      ("g:Profiler", "https://biit.cs.ut.ee/gprofiler/"),
                      ("Enrichr", "https://maayanlab.cloud/Enrichr/")):
        try:
            r = requests.get(url, timeout=30)
            print("service  reachable %-12s %s (HTTP %d)" % (name, url, r.status_code))
        except Exception as ex:
            print("service  NO ANSWER %-12s %s (%s)" % (name, url, str(ex)[:80])); ok = False
except ImportError:
    pass
sys.exit(0 if ok else 1)
PY
}

echo "repository : ${NANOAI_REPO}"
echo "data       : ${DATA}"
echo "results    : ${OUT}"
echo "python     : $(command -v "${PYTHON}" || echo "${PYTHON} (not found)")"

step "preflight"   preflight
step "ms_analysis" "${PYTHON}" "${HERE}/ms_analysis.py"

# ---------------------------------------------------------------- summary
echo
echo "================ ${SECTION}: summary ================"
printf '%s\n' "${SUMMARY[@]}"
echo "results: ${OUT}"
echo "logs   : ${LOGS}"
if [ "${NFAIL}" -gt 0 ]; then
  echo "${NFAIL} step(s) FAILED"
  echo "ms_analysis.py needs internet access (UniProt, g:Profiler, Enrichr) and the packages listed in README.md."
  exit 1
fi
echo "all steps PASSED"
