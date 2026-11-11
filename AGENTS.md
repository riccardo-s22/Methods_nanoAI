# AGENTS.md

Operating manual for a coding agent working in this repository. Every command below was run
from the repository root when this manual was written, with real values in place of `$ENVDIR`,
`REGEX` and `/path/to/...`. The exceptions are marked "not run".

## What this repository is

The public code and source data of a paper under review (Sirtori, Monroy Lopez, Roxbury,
Fallini): a Geneformer in-silico perturbation screen of a public human TBI snRNA-seq cohort,
linked to SWCNT nanosensor spectra of blast-injured human cortical organoids. Four section
directories (`01_` to `04_`) hold the original analysis scripts, `data/` holds their inputs, and
`reproduce.sh` reruns everything that can run on a CPU machine. The scripts are the as-run
analysis code with only their paths made portable; the as-run version is commit `2f6dd48`.

Figure numbers used here are the repository's: section 02 / "Fig. 2" is the Geneformer figure
(manuscript Figure 3), section 03 / "Fig. 3" is the protein corona figure (manuscript Figure 2).

## Setup

Python 3.13 with the pinned packages of `envs/spectral.yml`. Either of:

```bash
mamba env create -y -p "$ENVDIR" -f envs/spectral.yml          # ENVDIR: any path outside the repository
```

```bash
python3 -m venv "$ENVDIR" && "$ENVDIR/bin/pip" install -r envs/requirements-spectral.txt
```

Then either activate the environment or pass the interpreter explicitly, which every runner
accepts:

```bash
PYTHON="$ENVDIR/bin/python" bash reproduce.sh --list
```

Without `PYTHON`, the runners use `python3` from `PATH`. `reproduce.sh` prints the interpreter
and its package versions first, and stops with exit code 2 if NumPy, SciPy, scikit-learn, pandas
or matplotlib cannot be imported. `conda env create -f envs/spectral.yml` creates the same
environment under the name `nanoai-spectral` (not run here).

## Rerun everything

```bash
bash reproduce.sh --list        # every step, its class and the reason for each expected SKIP; runs nothing
bash reproduce.sh               # sections 01 02 03 04 in order, then repro/check_results.py; 1 h 48 min here
bash reproduce.sh --skip-slow   # without the permutation-heavy steps; 15 min here
bash reproduce.sh 02 03         # a subset, always run in ascending order
bash reproduce.sh --offline 02 03   # no internet: section 02 skips 5 steps, section 03 is not run
```

It works from any working directory and needs no environment variable. It ends with a table
(section, steps passed / failed / skipped, wall time) and exits 0 only if no step failed and no
regenerated table differs from `data/`. Expected on a fresh clone:

| Section | PASS | FAIL | SKIP |
|---|---|---|---|
| 01_swcnt_blast_classification | 33 | 0 | 4 |
| 02_geneformer_in_silico_screen | 16 | 0 | 13 |
| 03_protein_corona_proteomics | 2 | 0 | 0 |
| 04_drug_response_prediction | 17 | 0 | 7 |

`SKIP` is not a failure. Each SKIP line states its reason, and `--list` shows them in advance.

## Rerun one section

```bash
bash 01_swcnt_blast_classification/run_section.sh --only 'REGEX'    # also: --skip-slow, --skip-svg
bash 02_geneformer_in_silico_screen/run_section.sh --only 'REGEX'   # also: --offline, --with-project-tree
bash 03_protein_corona_proteomics/run_section.sh                    # no options; needs internet
bash 04_drug_response_prediction/run_section.sh --skip-slow --quiet # also: --sequential, --with-refit
```

`REGEX` is matched against the step names printed by `--list`. Section 04 has no `--only`: run
its scripts directly (they need no environment variable). Section 01 scripts run by hand need
`EXP2_TAG`, `EXP2_DIR` and `MIX_INVERTED=1` set together, so prefer its runner.

## Regenerate one figure or panel

Outputs are relative to `results/<section>/`. Times are wall times on an 8-core machine; for the
slow rows they are the step times of the final full run. The runners take the interpreter from
`PYTHON`; the rows that call `python3` directly use the `python3` on `PATH`, so activate the
environment first or write `"$ENVDIR/bin/python"` in place of `python3`.

| Panel | Command | Output | Time |
|---|---|---|---|
| Fig. 1a, 1b | `bash 01_swcnt_blast_classification/run_section.sh --only 'exp2_analysis2_24h'` | `exp2_24h_sham_tbi.png`; `svg/exp2_24h/exp2_24h_sham_tbi__panel_*.svg` | 3 min |
| Fig. 1c | `bash 01_swcnt_blast_classification/run_section.sh --only 'exp2_mixes_24h' && bash 01_swcnt_blast_classification/run_section.sh --only 'exp2_mix_intermediate_24h'` | `exp2_24h_dose_drivers.png`, `exp2_24h_mix_intermediate.png`; `svg/exp2_24h/` | 13 s |
| Fig. 1d | `bash 01_swcnt_blast_classification/run_section.sh --only 'exp2_classifier_compare_24h'` | `exp2_24h_classifier_dose_compare.png`; `svg/exp2_24h/` | 7 min |
| Fig. 2a | no runnable script: `figures/mt_cls_umap_RECONSTRUCTED.py` needs the MT checkpoint, which does not exist | none | |
| Fig. 2b | `bash 02_geneformer_in_silico_screen/run_section.sh --only 'crossmodel' && bash 02_geneformer_in_silico_screen/run_section.sh --only 'make_fig2_panelB'` | `drug_repurpose/crossmodel_consensus.tsv`; `figures/Figure2b_crossmodel_consensus.{svg,png}` | 2 s |
| Fig. 3 | `bash 03_protein_corona_proteomics/run_section.sh` | `figures/`, `TBI_blast_proteomics_results.xlsx`, `summary.json` | 2 to 5 min |
| Fig. 4a | no script: experimental timeline | none | |
| Fig. 4b | `python3 04_drug_response_prediction/exp3_drug_axis.py` | `exp3_3H_drug_axis.png` | 4 s |
| Fig. 4b as SVG | `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 shared/make_svg.py 04_drug_response_prediction/exp3_drug_axis.py 04_drug_response_prediction/svg EXP3_TP=3H` | `svg/exp3_3H_drug_axis__{full,panel_a,panel_b}.svg` | 5 s |
| Fig. 4c | no script: representative western blots | none | |
| Fig. 4d-f | no script: GraphPad bar charts of day-7 densitometry | none | |
| Fig. 4g | `python3 04_drug_response_prediction/exp3_3H_organoid_score.py` (scores only; the scatter was assembled outside any script) | `exp3_3H_organoid_scores.csv` | 2 s |
| Supp. Fig. 1 | `bash 01_swcnt_blast_classification/run_section.sh --only 'exp2_analysis2_2h' && bash 01_swcnt_blast_classification/run_section.sh --only 'exp2_mix_intermediate_2h' && bash 01_swcnt_blast_classification/run_section.sh --only 'exp2_dose_pca_figs_2h'` | `exp2_2h_sham_tbi.png` (a, b), `exp2_2h_mix_intermediate.png` (c), `exp2_2h_decision_dose_response.png` (d); `svg/exp2_2h/` | 4 min |
| Supp. Fig. 2 | `bash 01_swcnt_blast_classification/run_section.sh --only 'make_supp_auc'` | `supp_fig2_auc/supp_auc_results.json`, `supp_auc_blast_vs_sham.{svg,pdf,png}` | 41 min |
| Supp. Fig. 3 | `EXP2_PEAKFEATS_DIR=/path/to/folder bash 01_swcnt_blast_classification/run_section.sh --only 'exp2_chirality_classify$'` (the folder must hold `exp2_{2h,24h}_peakfeats.npz`, which are not in `data/`; without the variable the step prints SKIP) | `supp_fig3_chirality/exp2_chirality.{png,svg}`, `exp2_chirality_results.txt` | 6 s |
| Supp. Fig. 4 | `bash 01_swcnt_blast_classification/run_section.sh --only 'make_suppfig4'` (the Gaussian version prints SKIP unless `EXP2_PEAKFEATS_DIR` is set as for Supp. Fig. 3) | `supp_fig4_dilution/suppfig4_dilution_readouts.{svg,png,tif}`; with the variable also `suppfig4_gaussian_readouts.{svg,png,tif}` | 6 s |
| Supp. Fig. 5a, 5b | `bash 02_geneformer_in_silico_screen/run_section.sh --with-project-tree --only 'umap'` (prints SKIP unless the stage-01 AnnData and the FT embedding cache are under `GF_WORK`) | `umap/` | under 1 s when skipped |
| Supp. Fig. 5, metric panels | `bash 02_geneformer_in_silico_screen/run_section.sh --only 'make_suppfig5_model_metrics'` | `figures/suppfig5_metrics__*.{svg,png}` | 2 s |
| Supp. Fig. 6 | `bash 02_geneformer_in_silico_screen/run_section.sh --only 'excit_pathway_heatmap_[FM]T'` | `pathway_heatmap/excit_pathway_enrichment_heatmap_{FT,MT}.png`, `excit_pathway_enrichment_matrix_{FT,MT}.tsv` | 5 s |
| Supp. Fig. 7 | `MPLBACKEND=Agg python3 04_drug_response_prediction/exp3_perm_power.py` | `exp3_3H_perm_power.png` | 32 to 47 min |
| Supp. Fig. 7 as SVG | `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 shared/make_svg.py 04_drug_response_prediction/exp3_perm_power.py 04_drug_response_prediction/svg EXP3_TP=3H` (runs the script again; also the `exp3_perm_power__svg` step of a full section 04 run) | `svg/exp3_3H_perm_power__{full,panel_a,panel_b}.svg` | 31 min |
| Supp. Fig. 8 | same two commands as Fig. 4b | `exp3_3H_drug_axis.png`; `svg/exp3_3H_drug_axis__panel_{a,b}.svg` | |
| Supp. Fig. 9 | no script: GraphPad bar charts from western-blot densitometry | none | |
| Supp. Fig. 10a | no script: GraphPad | none | |
| Supp. Fig. 10b-c | `EXP3_WB_XLSX=/path/to/corr_wb_swcnt.xlsx python3 04_drug_response_prediction/figures/plot_wb_swcnt_corr.py` (the workbook is not in `data/`; without it the script exits 1 with "missing input") | `figures/wb_swcnt_corr/wb_swcnt_corr__{full,panel_a_NfL,panel_b_NfH}.svg` | 4 s |

The statistics quoted in the paper are printed to `results/<section>/logs/<step>.log` by the same
commands (direct `python3` commands print them to the terminal instead);
`docs/EXPECTED_RESULTS.md` names the log line for each.

## What cannot be rerun locally: do not try

- **GPU / cluster stages of section 02** (`00_environment`, `01_dataset_preparation`,
  `02_model_finetuning`, the perturbation arrays and stage 05 in `03_in_silico_perturbation`):
  218 GPU array tasks, about 865 GPU-hours, Slurm, a 0.9 GB GEO download and more than 8 GB RAM
  for preparation. Their per-gene output tables are shipped in `data/`; start from those.
- **Anything MT**: the MT checkpoint, the tokenized MT corpus and the MT per-gene tables do not
  exist anywhere. The seven `*_RECONSTRUCTED` scripts were rewritten from the Methods, were not
  used for the published results, and cannot run end to end. Do not try to repair or complete
  them, and do not present their output as the paper's.
- **Steps whose input is not in `data/`** (they print SKIP): Supp. Fig. 3 and the Gaussian
  Supp. Fig. 4 (`exp2_{2h,24h}_peakfeats.npz`), Supp. Fig. 10b-c (`corr_wb_swcnt.xlsx`), the
  6 h read (`exp3_load.py` self-check, `exp3_analysis.py`, `exp3_build_cache.py`), and in
  section 02 everything that needs the 20%-detection table or the stage-01 AnnData. Do not
  reconstruct, approximate or download substitutes for these inputs; ask the owner.
- **Refit of the chirality peaks from raw spectra** (`exp3_peakfit.py`, `exp3_peakfit_ratios.py`,
  `exp3_build_cache.py`): crashes under NumPy 2.5 because `numpy.trapz` no longer exists. Do not
  patch it. It needs the older stack of `envs/refit.yml`, passed as `PYTHON=...`.
- **Western blots, densitometry and GraphPad panels** (Fig. 4c-f, Supp. Figs. 9, 10a), the
  Fig. 4a timeline and the Fig. 4g scatter: no code exists.
- **The Spectronaut search** behind section 03: `data/` holds its exported report.

## Section 02 on a cluster

```bash
bash 02_geneformer_in_silico_screen/run_section.sh --print-cluster-plan   # prints the commands in run order; runs nothing
python3 02_geneformer_in_silico_screen/gf_paths.py                        # prints every resolved location and marks what is missing
```

Configuration is by environment variables, all documented at the top of
`02_geneformer_in_silico_screen/gf_env.sh`. The ones that matter:

| Variable | Default | Meaning |
|---|---|---|
| `GF_WORK` | `work/Geneformer_TBI` | project tree of the FT pipeline (h5ad, tokenized data, models, results) |
| `GF_CONDA_ENV` / `GF_PYTHON` | unset / `python` | environment with Geneformer and PyTorch (`envs/geneformer.yml`, then the `00_environment` scripts) |
| `GF_PREP_PYTHON` | `$GF_PYTHON` | interpreter with the scanpy / harmonypy / scrublet stack |
| `GF_FT_MODEL` | unset | a fine-tuned FT checkpoint to use instead of the pointer file |
| `GF_GEO` | `work/GSE209552` | the GEO download |

Slurm scripts are submitted from the repository root after creating
`results/02_geneformer_in_silico_screen/logs/`, because `#SBATCH --output` cannot expand
variables. They carry the partitions, constraints and limits of the original cluster: adapt
those lines only. They have no `--mail-user` line; pass one to `sbatch` if wanted. Once the tree
exists, `bash 02_geneformer_in_silico_screen/run_section.sh --with-project-tree` runs the CPU
steps that read it. None of the cluster commands was run for this manual.

## Where outputs and logs land

| What | Where |
|---|---|
| inputs, read-only | `data/<section>/` (override: `NANOAI_DATA`) |
| outputs | `results/<section>/` (override: `NANOAI_RESULTS`) |
| one log per step | `results/<section>/logs/<step>.log` |
| section 04 status table | `results/04_drug_response_prediction/logs/run_status.tsv` |
| section 02 stage 07 and 09 tables | `work/Geneformer_TBI/results/` (override: `GF_WORK`) |
| `reproduce.sh` | `results/reproduce/<section>.log`, `summary.tsv`, `check_results.txt` |

`results/` and `work/` are git-ignored. Nothing is written to the current directory.

## Verify a rerun, or explain a number that differs

```bash
python3 repro/check_results.py                         # regenerated tables against data/: PASS / FAIL / NOT PRODUCED
python3 04_drug_response_prediction/verify_outputs.py  # section 04: tables and printed statistics against the published run
```

1. Find the statistic in `docs/EXPECTED_RESULTS.md`. It gives the paper value, the value a rerun
   gives, and the log line or file column. Known mismatches with the manuscript are listed per
   section under **Differences**: if yours is there, it is expected, and not yours to fix.
2. Compare the interpreter line at the top of the `reproduce.sh` output with `envs/spectral.yml`.
3. Check that you are reading the right variant: 2 h or 24 h read; floored
   (`results/04_.../`) or true-log (`results/04_.../sensitivity_true_log_ratios/`) features.
4. Rerun the single step and read its log.
5. Report both values, the script and the log line. Do not change anything to close the gap.

## Known caveats

- **Well A1 is a placeholder.** In the section 01 spectra, column 1 (well A1) is a synthetic
  0 to 1000 ramp, not data; `exp2_load.py` excludes it, leaving 95 wells. The section 04 files
  have a similar ramp column, excluded by `exp3_load.py`.
- **`MIX_INVERTED=1` is the published dilution labelling** (column 10 = 25% blast, column 12 =
  75%). The runner and the scripts that use the dilution wells set it; `exp2_load.py` alone does not.
- **`EXP2_TAG` is only a label.** `EXP2_DIR` selects the spectra, and its default is the 2 h read.
  `EXP2_TAG=24h` without `EXP2_DIR` analyses 2 h spectra under a 24 h title.
- **Corrected blast-plate layout (section 04).** On the blast plate, columns 1-4 are MG132,
  5-8 vehicle, 9-12 resveratrol. `exp3_load.py` encodes this; `README_plate_layout_original.txt`
  has resveratrol and MG132 swapped and must not be used. For the same reason the two block
  titles printed by `supporting_checks/exp3_drug_swap.py` are exchanged.
- **The cached fits are authoritative.** Section 04 reads
  `data/04_drug_response_prediction/chirality_features/exp3_3H_peakfeats.npz`. Refitting crashes
  under the Python 3.13 stack of `envs/spectral.yml`, and with `numpy.trapz` restored it returns
  different features (ratio AUC 0.934, not 0.938). The fits were made with the Python 3.10 stack
  of `envs/refit.yml`, which regenerates the cache bit for bit.
- **Floored vs true log ratios.** Every published chirality number uses
  `log10((H_i + 1e-6) / (H_j + 1e-6))`, although the Methods define `log10(H_i / H_j)`.
  `sensitivity_true_log_ratios/` reruns with the true ratio (AUC 0.887 instead of 0.938). Both
  write a file named `exp3_3H_organoid_scores.csv`, in different folders.
- **Coverage gate is not the eligibility rule.** `*_filtered_high_confidence_gene_ranked.tsv`
  (12,689 to 13,393 genes) is a coverage gate. The paper's eligible universe is 7,374 genes
  (detected in at least 20% of either group), applied downstream, and it needs a detection
  table that is not in `data/`.
- **Eight rows, six consensus drugs.** `crossmodel_consensus.tsv` has eight rows; `vx` and
  `sodium` are fragments left by the name parser and are dropped by the figure scripts.
- **Pass-through in section 02.** In a default run the two Supp. Fig. 6 matrices are read from
  `data/` and written back, so their `SAME` lines are not an independent recomputation. They
  hold FDR-corrected p values (the heatmap scripts take -log10 when drawing).
- **Section 03 depends on live services.** Gene symbols (5 of 6,246 so far) and the
  `3_GSEA_preranked` sheet drift with UniProt and Enrichr releases; `summary.json` and every
  reported number should stay identical. Label positions in the volcano plot vary between runs.
- **Unseeded parts, left as published:** `gap` values of the elastic-net models in
  `supp_auc_results.json`, the Supp. Fig. 4 dilution image, the x jitter of the strip plots in
  Fig. 4b. Images differ between runs; printed statistics do not.
- **`shared/make_svg.py` needs the thread variables in the shell**
  (`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`), as in the Fig. 4b row. Passed
  as `ENV=VAL` arguments they come too late and the run is many times slower.
- **Arial.** The figure scripts ask for Arial and fall back to DejaVu Sans when it is missing,
  which changes text metrics, not content.
- **`data/.../chirality_features/exp2_*_mix_peakfeats.npz`** hold the dilution and
  negative-control wells only; the sham/blast fits are the files that are not shipped.

## Rules

1. Never edit analysis logic to "make it run": no statistic, model, hyperparameter, seed,
   permutation count, filter, threshold, feature definition, label or plot content. Allowed in
   existing scripts: where inputs are read, where outputs are written, and the `sys.path` bootstrap.
2. Never write, move, rename, add or delete anything under `data/`. Outputs go under `results/`
   (or `work/` for section 02).
3. When a result differs from the paper or from a shipped table, report it with both values.
   Never tune a parameter, seed, tolerance or package version toward the published value.
4. Do not commit `results/` or `work/`. Do not weaken `.gitignore`.
5. Do not shorten a run by lowering permutation counts in the repository; use `--skip-slow`
   and say which steps were skipped.
6. Do not add data files, substitute inputs or mock services to turn a SKIP into a PASS.
7. Keep machine-specific paths, user names and e-mail addresses out of the repository; use the
   environment variables instead.
8. If a step fails, read `results/<section>/logs/<step>.log` before changing anything, and fix
   the environment or the invocation, not the analysis.
