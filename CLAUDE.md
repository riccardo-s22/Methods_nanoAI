# CLAUDE.md

@AGENTS.md

## Specific to Claude Code

- Start with `bash reproduce.sh --list`: it runs nothing and shows which steps run locally and
  which will print `SKIP`, with the reason.
- `bash reproduce.sh` takes about 1 h 50 min. Start it with `run_in_background` and do not poll
  it in a loop; read `results/reproduce/summary.tsv` and `results/reproduce/check_results.txt`
  when it has exited. The same applies to `01_swcnt_blast_classification/run_section.sh`
  and `04_drug_response_prediction/run_section.sh` without `--skip-slow`.
- Sections 02 and 03 query web services (Enrichr, STRING, UniProt, g:Profiler). In a sandbox
  without network access use `bash reproduce.sh --offline` and say that section 03 was not run.
- Do not add permission rules or hooks that allow writing under `data/`.
- When a number differs from the paper, look it up in `docs/EXPECTED_RESULTS.md` first: the known
  differences are listed there with both values.
