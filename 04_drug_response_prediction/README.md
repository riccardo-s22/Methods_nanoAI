# Results 4 — early prediction of a drug's neuroprotective effect

**Figure 4; Supplementary Figures 7–10.**

Organoids were blast- or sham-exposed, treated for 12 h with 100 µM resveratrol
(nominated by the Geneformer screen), 10 µM MG132 or vehicle, and their
conditioned medium read with the SWCNT sensor 36 h after injury. The same
organoids were then collected on day 7 for western blotting, so the optical
readout can be tested against the biochemical outcome it is meant to anticipate.

## Rerun in one command

    bash 04_drug_response_prediction/run_section.sh

It works from any working directory, needs no environment variable, reads only
`data/04_drug_response_prediction/` and writes only under
`results/04_drug_response_prediction/` (git-ignored). It runs every step below
that can be rerun from the shipped data, in README order (the four slow steps
are started first and run in the background, see the table), tees each step's
output to `results/04_drug_response_prediction/logs/<step>.log`, prints one
`PASS`, `FAIL` or `SKIP` line per step, finishes with `verify_outputs.py`
(regenerated tables and printed statistics against the published ones) and
exits non-zero if any step failed. A machine-readable status table is written
to `logs/run_status.tsv`.

| Option | Effect |
|---|---|
| *(none)* | everything that can be reproduced from the shipped data: 17 steps run, 7 are skipped for the reasons given below. The four slow steps (`exp3_perm_power.py`, its true-log copy, and the SVG re-render of each; 20,000 permutations per feature set) are started first and run concurrently in the background while the others proceed. The two complete runs timed here took 63 min and 92 min, both on a machine busy with other jobs; with 4 idle cores it should come close to the runtime of one `exp3_perm_power.py` run (32 min when run alone). See *Runtimes*. |
| `--skip-slow` | skips those four slow steps. About 10 min. |
| `--sequential` | runs the slow steps one after the other, in README order. This is also what happens on a machine with fewer than 4 CPUs; `--parallel` forces the concurrent mode. |
| `--quiet` | script output goes to the logs only; the terminal shows the status lines. |
| `--with-refit` | also attempts the three scripts that refit the peaks from the raw spectra. They crash under the NumPy version of the paper; see *Known limitations*. |

The same switches exist as environment variables (`SKIP_SLOW=1`, `PARALLEL=0`
or `1`, `QUIET=1`, `WITH_REFIT=1`). `PYTHON` selects the interpreter (default `python3`).

Software: Python 3.13, NumPy 2.5, SciPy 1.18, scikit-learn 1.9, pandas 3.0,
matplotlib 3.11; `openpyxl` only for `plot_wb_swcnt_corr.py`.

Steps that are **not** run from a fresh clone, each printed as one `SKIP` line:

| Step | Why | How to enable |
|---|---|---|
| `plot_wb_swcnt_corr.py`, `plot_wb_swcnt_corr_truelog.py` (Supp. Fig. 10b–c) | the western-blot densitometry workbook `corr_wb_swcnt.xlsx` is not in the repository (available from the corresponding author) | `EXP3_WB_XLSX=/path/to/corr_wb_swcnt.xlsx`, or place it at `data/04_drug_response_prediction/western_blot/corr_wb_swcnt.xlsx` |
| `exp3_analysis.py`, `exp3_load.py` run as a script | both loop over the 3 h **and** the 6 h read; the 6 h spectra are not in the repository | `EXP3_DIR=/folder/with/both/reads` |
| `exp3_peakfit.py`, `exp3_peakfit_ratios.py`, `exp3_build_cache.py` | refit from raw spectra; crashes under NumPy 2.5 (see *Known limitations*). `exp3_build_cache.py` also needs the 6 h read | run them with an interpreter that has SciPy 1.15 and a NumPy that still has `numpy.trapz` (see *Known limitations*): `PYTHON=/path/to/that/python bash 04_drug_response_prediction/run_section.sh` runs them automatically. `--with-refit` attempts them with the current interpreter regardless |

## Design

Two 96-well plates, one sham and one blast, read 3 h after SWCNT addition.

| Plate | Columns 1–4 | Columns 5–8 | Columns 9–12 |
|---|---|---|---|
| Sham | sham vehicle | resveratrol | MG132 |
| Blast | MG132 | blast vehicle | resveratrol |

32 wells per condition, except MG132-on-sham, which is 24 (rows A–F only;
`README_plate_layout_original.txt` records the exceptions — wells G9–12 and
H9–10 are SWCNT in unconditioned medium, H11–12 medium without SWCNTs).
184 wells in total.

**The blast plate's resveratrol and MG132 columns are swapped relative to the
original plate README.** This was confirmed against lab records, and
`exp3_load.py` encodes the corrected assignment. `supporting_checks/exp3_drug_swap.py`
is the script that tested both labellings side by side so the correct one could
be chosen from records rather than from which fit the hypothesis.

Each well is reduced to its fitted emission peaks: a split Gaussian per
canonical (n,m) species, giving 12 peak-area fractions and 66 pairwise log10
intensity ratios.

## Run order

| Script | What it does |
|---|---|
| `exp3_load.py` | Loader. Five excitations × 512 emission rows × 98 columns per plate and time point; assigns the corrected plate layout. |
| `exp3_peakfit.py` | Split-Gaussian fits of the 12 canonical chiralities per well; the 12 peak-area fractions. |
| `exp3_peakfit_ratios.py` | Adds the 66 scale-invariant pairwise log10 intensity ratios. |
| `exp3_build_cache.py` | Rebuilds the `exp3_{TP}_peakfeats.npz` feature caches. Run after any plate-label change. |
| `exp3_export_peakfeats.py` | Exports the cached features to CSV (per well: group, plate row/column, 12 `gauss_max`, 12 `gauss_auc`, 12 area fractions, 66 log ratios). |
| `exp3_analysis.py` | Sham-versus-blast discrimination and the drug arms; PERMANOVA. |
| `exp3_drug_axis.py` | **Figure 4b** and **Supplementary Figure 8**. Fits the LDA sham→blast disease axis on vehicle-treated wells only, scores every other well without refitting, and decomposes each drug effect vector into its along-axis and orthogonal components. |
| `exp3_drug_effects.py` | The 2 × 3 factorial (condition × drug), including the difference-in-differences contrast for drug effect specific to blast exposure. |
| `exp3_3H_organoid_score.py` | **Figure 4g** (the scores; the scatter against day-7 NfL was assembled outside the script) and the input to `figures/plot_wb_swcnt_corr.py`. Per-organoid injury score for correlation with day-7 biochemistry. Sham and blast wells are scored leave-one-out (standardization and LDA refit without the scored well); drug-treated wells use the model fitted on all sham and blast wells. Scores are rescaled to `fc = 1 + (s − μ_sham)/(μ_blast − μ_sham)` so sham = 1 and blast = 2. |
| `exp3_perm_power.py` | **Supplementary Figure 7**. Permutation nulls for both feature sets — the 2,560-dimensional shape vectors reduced to 20 PCs, and the 66 chirality log ratios — and how many permutations are needed to certify significance. |

### Commands, inputs, outputs and runtimes

No environment is needed: every script locates the repository from its own
path, so it runs from any working directory with no variable set. Commands
below are written from the repository root. `D` stands for
`data/04_drug_response_prediction` and `R` for `results/04_drug_response_prediction`.

| Script | Command | Reads | Writes | Runtime |
|---|---|---|---|---|
| `exp3_load.py` | `python3 04_drug_response_prediction/exp3_load.py` | `D/spectra/*.txt`, plus the 6 h read (not shipped) | nothing; prints a per-group summary of each plate and read | 1 s with both reads. With the shipped data it prints the sham-plate 3 h block and stops with `AssertionError: missing excitations for SHAM 6H`; the runner skips it. As a module it is imported by every other script and needs only the 3 h files. |
| `exp3_peakfit.py` | `python3 04_drug_response_prediction/exp3_peakfit.py` | `D/spectra/*.txt` | `R/exp3_3H_peakfit.png` | crashes after about 2 s under NumPy 2.5 (see *Known limitations*); skipped by the runner |
| `exp3_peakfit_ratios.py` | `python3 04_drug_response_prediction/exp3_peakfit_ratios.py` | `D/spectra/*.txt` | `R/exp3_3H_peakfit_ratios.png` | same |
| `exp3_build_cache.py` | `python3 04_drug_response_prediction/exp3_build_cache.py` | `D/spectra/*.txt`, plus the 6 h read (not shipped) | `R/rebuilt_cache/exp3_{3H,6H}_peakfeats.npz` — never the shipped cache in `D/chirality_features/` | same |
| `exp3_export_peakfeats.py` | `python3 04_drug_response_prediction/exp3_export_peakfeats.py` | `D/chirality_features/exp3_3H_peakfeats.npz` | `R/exp3_3H_chirality_peak_features.csv` (190 wells × 107 columns). The line `missing …exp3_6H_peakfeats.npz` is expected: no 6 h cache is shipped. | under 1 s |
| `exp3_analysis.py` | `EXP3_DIR=/folder/with/both/reads python3 04_drug_response_prediction/exp3_analysis.py` | 3 h and 6 h spectra (6 h not shipped) | `R/exp3_sham_tbi_drugs.png` | 3.5 min with both reads; skipped by the runner otherwise |
| `exp3_drug_axis.py` | `python3 04_drug_response_prediction/exp3_drug_axis.py` | `D/chirality_features/exp3_3H_peakfeats.npz` | `R/exp3_3H_drug_axis.png` | 6 s |
| `exp3_drug_axis.py`, as SVG panels | `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 shared/make_svg.py 04_drug_response_prediction/exp3_drug_axis.py 04_drug_response_prediction/svg EXP3_TP=3H` | same | `R/svg/exp3_3H_drug_axis__{full,panel_a,panel_b}.svg` (panel a = Fig. 4b / Supp. Fig. 8A, panel b = Supp. Fig. 8B) | 4 s |
| `exp3_drug_effects.py` | `python3 04_drug_response_prediction/exp3_drug_effects.py` | `D/chirality_features/exp3_3H_peakfeats.npz` | `R/exp3_3H_drug_effects.png` | 3 s |
| `exp3_3H_organoid_score.py` | `python3 04_drug_response_prediction/exp3_3H_organoid_score.py` | `D/chirality_features/exp3_3H_peakfeats.npz` | `R/exp3_3H_organoid_scores.csv` (reference copy: `D/per_well_scores/exp3_3H_organoid_scores.csv`) | 4 s |
| `exp3_perm_power.py` | `python3 04_drug_response_prediction/exp3_perm_power.py` | `D/spectra/*.txt`, `D/chirality_features/exp3_3H_peakfeats.npz` | `R/exp3_3H_perm_power.png` | 32 min — **slow** |
| `exp3_perm_power.py`, as SVG panels | `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 shared/make_svg.py 04_drug_response_prediction/exp3_perm_power.py 04_drug_response_prediction/svg EXP3_TP=3H` | same | `R/svg/exp3_3H_perm_power__{full,panel_a,panel_b}.svg` (Supp. Fig. 7A, 7B) | 32 min — **slow** (the script is run a second time) |
| `verify_outputs.py` | `python3 04_drug_response_prediction/verify_outputs.py` | `R/*.csv`, `R/logs/*.log`, the reference tables in `D/` | nothing; prints `OK`, `MISMATCH` or `NOT RUN` per check | 1 s |

`shared/make_svg.py` re-runs a figure script unchanged and intercepts
`savefig`, so the SVG steps print the same statistics as the plain runs but
write SVG instead of PNG. Its second argument is the output folder, relative to
the results tree. The published panels of Fig. 4b and Supp. Figs. 7 and 8 are
these SVGs. The `saved ….png` line in the log of an SVG run is the script's own
message; that run writes the SVG files only. Set the three thread variables in
the shell for a `make_svg.py` call, as in the commands above and as the runner
does: the scripts pin BLAS to one thread themselves, but under the harness NumPy
is already loaded by then (see [`../shared/README.md`](../shared/README.md)).

The dots of the strip plots (`exp3_drug_axis.py`, `exp3_drug_effects.py`,
`exp3_analysis.py`, `exp3_peakfit*.py`) are spread horizontally with an
unseeded random jitter, so their x positions differ from run to run. Their y
values, the group means and every printed statistic do not.

All scripts read the feature cache from `D/chirality_features/`; none writes
there. Outputs of a script can be redirected with `EXP3_OUT` (see
*Environment variables*).

## figures/

`plot_wb_swcnt_corr.py` — **Supplementary Figure 10b–c**. Correlates the mean SWCNT score of each of the six conditions (sham or blast ×
vehicle, resveratrol, MG132) against the mean NfL and NfH western-blot values of
the same condition, normalized to GAPDH and expressed relative to vehicle-treated
sham, by Pearson and Spearman (n = 6 conditions). Data: `corr_wb_swcnt.xlsx`.

Figure 4a (timeline), 4c (representative blots) and 4d–f (day-7 p-tau, NfH and
NfL quantification), together with Supplementary Figures 9 and 10a, involve no
code: blots were quantified in ImageJ and the bar charts drawn in GraphPad Prism
(two-way ANOVA with Tukey's multiple-comparisons test).

| Script | Command | Reads | Writes | Runtime |
|---|---|---|---|---|
| `figures/plot_wb_swcnt_corr.py` | `EXP3_WB_XLSX=/path/to/corr_wb_swcnt.xlsx python3 04_drug_response_prediction/figures/plot_wb_swcnt_corr.py` | `corr_wb_swcnt.xlsx`, sheet `Sheet2` (**not in the repository**); the per-well scores `R/exp3_3H_organoid_scores.csv` if `exp3_3H_organoid_score.py` has been run, otherwise the shipped `D/per_well_scores/exp3_3H_organoid_scores.csv` | `R/figures/wb_swcnt_corr/wb_swcnt_corr__{full,panel_a_NfL,panel_b_NfH}.svg` | 4 s; skipped by the runner when the workbook is absent |

Without the workbook the script stops with a `missing input:` message and exit
status 1. Its default location is
`data/04_drug_response_prediction/western_blot/corr_wb_swcnt.xlsx`; the script
prints which workbook and which score table it read.

## supporting_checks/

| Script | Question |
|---|---|
| `exp3_3H_pressure.py` | Does the headline sham-versus-blast result survive the spatial plate gradient and outlier wells? Unlike the section 1 plate, sham and blast here sit on different plates, so position is not balanced by design and has to be checked. |
| `exp3_drug_swap.py` | Tests both candidate labellings of the blast plate's drug columns. The sham→blast axis is unchanged by the swap, since it is built from sham and blast wells, neither of which is affected. |

| Script | Command | Reads | Writes | Runtime |
|---|---|---|---|---|
| `supporting_checks/exp3_3H_pressure.py` | `python3 04_drug_response_prediction/supporting_checks/exp3_3H_pressure.py` | `D/spectra/*.txt` | nothing; prints | 33 s |
| `supporting_checks/exp3_drug_swap.py` | `python3 04_drug_response_prediction/supporting_checks/exp3_drug_swap.py` | `D/chirality_features/exp3_3H_peakfeats.npz` | nothing; prints | 2 s |

Two things to know when reading their output:

- `exp3_3H_pressure.py` part [E] refits the peaks from the raw spectra (it
  restores `numpy.trapz` itself, so it runs) and prints a 66-ratio AUC of
  **0.934**, against the published 0.938 obtained from the cache. The same 0.934
  is in the original log of this script. See *Known limitations* on why a refit
  does not return the cached values exactly.
- `exp3_drug_swap.py` was written before the cache was rebuilt with the
  corrected labels, and its printed block titles and column annotations
  describe that earlier cache. Run on the shipped (corrected) cache, the block
  titled `ORIGINAL (README)` therefore holds the **corrected** assignment
  (resveratrol = blast-plate columns 9–12; 43% rescue) and the block titled
  `SWAPPED` holds the assignment of the original plate README (5%). The numbers
  are the same as in the original run; only the two titles are exchanged.

## sensitivity_true_log_ratios/

Re-runs the chirality analyses with true `log10(H_i/H_j)` ratios instead of
`log10((H_i + 1e-6)/(H_j + 1e-6))` — see the repository README for why the floor
matters. `README_truelog.md` tabulates old against new for every affected
statistic; `run_all.sh` drives the re-run.

This section is where the floor actually bites:

| Statistic | Floored (published) | True log ratio |
|---|---|---|
| LOO AUC, 66 ratios | 0.938, z = 4.5 | 0.887, z = 4.0 |
| LOO AUC, full shape | 0.903, z = 4.3 | 0.903, z = 4.3 (unaffected) |

`exp3_peakfit_ratios_truelog.py` reproduces the published floored values exactly
when run with `FLO=1e-6`, which is what licenses the comparison. It reads the
cached peak fits rather than refitting from raw spectra, because refitting fails
in the current SciPy — see *Known limitations*.

    bash 04_drug_response_prediction/sensitivity_true_log_ratios/run_all.sh

runs the re-run on its own (options `--skip-slow`, `--sequential`, `--quiet`,
as for `run_section.sh`, which also calls it). `RT` stands for `R/sensitivity_true_log_ratios`.

| Script | Command | Reads | Writes | Runtime |
|---|---|---|---|---|
| `exp3_3H_organoid_score_truelog.py` | `python3 04_drug_response_prediction/sensitivity_true_log_ratios/exp3_3H_organoid_score_truelog.py` | `D/chirality_features/exp3_3H_peakfeats.npz` | `RT/exp3_3H_organoid_scores.csv` (reference copy: `D/per_well_scores/exp3_3H_organoid_scores_truelog.csv`) | 1 s |
| `exp3_drug_axis_truelog.py` | `python3 04_drug_response_prediction/sensitivity_true_log_ratios/exp3_drug_axis_truelog.py` | same cache | `RT/exp3_3H_drug_axis.png`; as SVG panels `RT/svg/exp3_3H_drug_axis__*.svg` (`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 shared/make_svg.py 04_drug_response_prediction/sensitivity_true_log_ratios/exp3_drug_axis_truelog.py 04_drug_response_prediction/sensitivity_true_log_ratios/svg EXP3_TP=3H`) | 4 s, and 4 s for the SVG run |
| `exp3_peakfit_ratios_truelog.py` | `python3 04_drug_response_prediction/sensitivity_true_log_ratios/exp3_peakfit_ratios_truelog.py` | same cache, `D/spectra/*.txt` | `RT/exp3_3H_peakfit_ratios.png` | 4.6 min |
| the same with the published floor | `FLO=1e-6 EXP3_OUT=results/04_drug_response_prediction/sensitivity_true_log_ratios/checkfloor_FLO1e-6 python3 04_drug_response_prediction/sensitivity_true_log_ratios/exp3_peakfit_ratios_truelog.py` | same | `RT/checkfloor_FLO1e-6/exp3_3H_peakfit_ratios.png`; log `R/logs/exp3_peakfit_ratios_truelog__FLO1e-6.log` | 4.5 min |
| `exp3_perm_power_truelog.py` | `python3 04_drug_response_prediction/sensitivity_true_log_ratios/exp3_perm_power_truelog.py` | same cache, `D/spectra/*.txt` | `RT/exp3_3H_perm_power.png`; as SVG panels `RT/svg/exp3_3H_perm_power__*.svg` (same `make_svg.py` form) | 32 min, and 32 min for the SVG run — **slow** |
| `plot_wb_swcnt_corr_truelog.py` | `EXP3_WB_XLSX=/path/to/corr_wb_swcnt.xlsx python3 04_drug_response_prediction/sensitivity_true_log_ratios/plot_wb_swcnt_corr_truelog.py` | `corr_wb_swcnt.xlsx` (**not in the repository**); `RT/exp3_3H_organoid_scores.csv` if present, otherwise `D/per_well_scores/exp3_3H_organoid_scores_truelog.csv` | `RT/svg/wb_swcnt_corr__{full,panel_a_NfL,panel_b_NfH}.svg` | 3 s; skipped by the runner when the workbook is absent |

The true-log scripts are copies of the published ones in which only the floor
(`flo = 0.0`) and the output folder differ. Their output file names are the
same as those of the published scripts, so they are kept apart by folder:
everything from this directory lands in `RT/`.

## Inputs

| File | What it is |
|---|---|
| `SHAM_EX_3_3H__{570,640,670,750,780}.txt`, blast equivalents | Raw NIR emission, one file per excitation per plate per time point. |
| `exp3_3H_peakfeats.npz` | Cached split-Gaussian fits with corrected labels. |
| `corr_wb_swcnt.xlsx` | Western-blot densitometry for the SWCNT-versus-biochemistry correlation. |

Where these are:

| Input | Location | Override |
|---|---|---|
| 3 h spectra, 10 files | `data/04_drug_response_prediction/spectra/` | `EXP3_DIR` |
| `exp3_3H_peakfeats.npz` | `data/04_drug_response_prediction/chirality_features/` | `EXP3_CACHE_DIR` |
| published per-well scores (reference copies) | `data/04_drug_response_prediction/per_well_scores/` | — |
| `corr_wb_swcnt.xlsx` | **not in the repository**; expected at `data/04_drug_response_prediction/western_blot/corr_wb_swcnt.xlsx` | `EXP3_WB_XLSX` |
| 6 h spectra (`SHAM_EX3_6H__*.txt`, `TBI_EX3_6H__*.txt`) | **not in the repository**; read from the same folder as the 3 h files | `EXP3_DIR` |

`data/04_drug_response_prediction/per_well_scores/exp3_3H_organoid_scores_shape.csv`
has no generating script in this directory.

## Environment variables

None is required. All of them only say where things are; no analysis parameter
is read from the environment except `EXP3_TP` and `FLO`, which the original
scripts already read.

| Variable | Meaning | Default |
|---|---|---|
| `NANOAI_DATA` | source-data tree | `<repo>/data` |
| `NANOAI_RESULTS` | results tree | `<repo>/results` |
| `EXP3_DIR` | folder holding the spectra `.txt` files | `$NANOAI_DATA/04_drug_response_prediction/spectra` |
| `EXP3_TP` | read analysed by `exp3_drug_axis.py` and `exp3_drug_effects.py` | `3H` (the only read whose cache is shipped) |
| `EXP3_CACHE_DIR` | folder holding `exp3_{TP}_peakfeats.npz` | `$NANOAI_DATA/04_drug_response_prediction/chirality_features` |
| `EXP3_OUT` | output folder of the script being run | `$NANOAI_RESULTS/04_drug_response_prediction` (true-log scripts: its `sensitivity_true_log_ratios/` subfolder; `exp3_build_cache.py`: `rebuilt_cache/`; `plot_wb_swcnt_corr.py`: `figures/wb_swcnt_corr/`) |
| `EXP3_WB_XLSX` | western-blot densitometry workbook | `$NANOAI_DATA/04_drug_response_prediction/western_blot/corr_wb_swcnt.xlsx` |
| `EXP3_SCORES_CSV` | per-well score table read by the two `plot_wb_swcnt_corr*` scripts | the regenerated table under `results/` if present, otherwise the shipped one |
| `FLO` | floor added to the peak intensities in `exp3_peakfit_ratios_truelog.py` | `0.0` (`1e-6` gives back the published values) |

## vendored/

Two modules the scripts import that were not in the repository before; see
`vendored/README.md` for where each comes from.

| File | Used by | For |
|---|---|---|
| `vendored/chirality_interpolated_extraction.py` | `exp3_peakfit.py`, `exp3_peakfit_ratios.py`, `exp3_build_cache.py`, `supporting_checks/exp3_3H_pressure.py`, `sensitivity_true_log_ratios/exp3_peakfit_ratios_truelog.py` | the split-Gaussian peak extractor and the 12 reference (n,m) coordinates |
| `vendored/exp3_7d_load.py` | `supporting_checks/exp3_3H_pressure.py` | `outlier_keep()` only |

## Known limitations

- **The peak fitting is not currently reproducible from raw spectra.** The
  pipeline loads the cached `exp3_3H_peakfeats.npz`; refitting the split
  Gaussians from the raw `.txt` files fails under the installed SciPy. The cache
  is the authoritative feature matrix for this section.
- `gauss_auc` falls back to trapezoidal integration when a Gaussian fit fails,
  so a small number of area values are not fitted areas.
- The 66 pairwise log ratios of 12 quantities have rank 11, so a PCA asking for
  more than 11 components returns noise directions.

What happens exactly when the refit is attempted (Python 3.13.13, NumPy 2.5.1,
SciPy 1.18.0):

- `exp3_peakfit.py`, `exp3_peakfit_ratios.py` and `exp3_build_cache.py` stop
  after about 2 s with `AttributeError: module 'numpy' has no attribute 'trapz'`.
  The sequence is: for some wells SciPy's `curve_fit` raises
  `ValueError: Initial guess is outside of provided bounds`; the extractor
  catches that and takes its documented fallback (interpolated maximum as peak
  height, trapezoidal area); the fallback calls `numpy.trapz`, which NumPy 2.5
  no longer has. So it is the fallback, not the fit, that brings the script
  down. The scripts are left as they are.
- `supporting_checks/exp3_3H_pressure.py` defines `numpy.trapz` itself before
  refitting, so its refit runs. It gives a 66-ratio AUC of 0.934 rather than
  the cached 0.938.
- A refit in this environment with `numpy.trapz` restored does **not** return
  the shipped cache exactly. Over the 190 wells × 12 species of the 3 h read,
  883 of 2,280 peak heights are bit-identical, 1,824 agree to a relative 1e-4,
  and the largest relative difference is 1.1%; for the areas (`gauss_auc`),
  336 of 2,280 differ by more than 10%.
- With **Python 3.10.20, NumPy 2.2.6, SciPy 1.15.3 and scikit-learn 1.7.2** the
  same unmodified scripts do run from the raw spectra, and they return the
  published values exactly: `exp3_build_cache.py` regenerates the shipped
  `exp3_3H_peakfeats.npz` bit for bit (all 4,560 values; it needs the 6 h read
  to run to completion), and `exp3_peakfit_ratios.py` prints 0.613 / 0.938 /
  0.920 with the published top ratios ((8,3)/(10,2) d = −0.86, (8,3)/(8,4)
  d = −0.76, (6,5)/(7,6) d = +0.77), and `exp3_peakfit.py` the 12-fraction AUC
  of 0.612 (permutation p = 0.0780) of its original run. So the cache is what the fitting code
  gives under that older SciPy and NumPy, and the peak fitting is reproducible
  from the raw spectra in such an environment; it is the newer stack that
  departs from it. Which change between the two stacks is responsible was not
  investigated.
- `exp3_load.py` run as a script and `exp3_analysis.py` need the 6 h read and
  stop at `AssertionError: missing excitations for SHAM 6H` without it.
- `exp3_3H_peakfeats.npz` holds 190 wells: the 184 analysed organoid wells and
  the 6 wells of SWCNT in unconditioned medium.

## Reproduction status

Rerun on the shipped data with `run_section.sh`; `verify_outputs.py` repeats
these comparisons.

| Output | Rerun against the shipped or published value |
|---|---|
| `exp3_3H_organoid_scores.csv` | identical to `data/…/per_well_scores/exp3_3H_organoid_scores.csv` (184 rows, same order and labels, max abs diff 0) |
| true-log `exp3_3H_organoid_scores.csv` | identical to `data/…/per_well_scores/exp3_3H_organoid_scores_truelog.csv` (max abs diff 0) |
| `exp3_3H_chirality_peak_features.csv` | identical to `data/…/chirality_features/exp3_3H_chirality_peak_features.csv` (190 × 107, max abs diff 0) |
| Supp. Fig. 7 (`exp3_perm_power.py`) | AUC 0.9033 (shape) and 0.9375 (66 ratios); null 95th percentile 0.636 and 0.635; 0 of 20,000 permutations reach the observed AUC; z = 4.32 and 4.54 — as published |
| Fig. 4b, Supp. Fig. 8 (`exp3_drug_axis.py`) | every printed line equals the log of the published run: scores −0.220 ± 0.111 (sham), 0.220 ± 0.124 (blast), 0.029 ± 0.187 (blast + resveratrol, p = 0.0001), 0.199 ± 0.133 (blast + MG132, p = 0.5267); rescue 43% (28–60%) and 5% (−10 to 18%); difference-in-differences −0.041 (p = 0.4796) and −0.148 (p = 0.0008); orthogonal components 4.03, 1.29, 4.33, 3.84 |
| Supp. Fig. 10b–c (`plot_wb_swcnt_corr.py`) | not rerun from the repository (workbook not shipped). With the workbook: NfL r = −0.809 (p = 0.0514), ρ = −0.829 (p = 0.0416); NfH r = −0.513 (p = 0.298), ρ = −0.543 (p = 0.266) — as published |
| true-log re-run | 66-ratio AUC 0.8867, z = 4.02; per-well score AUC 0.913; `FLO=1e-6` gives back 0.613 / 0.938 / 0.920 — as in `README_truelog.md` |
| Fig. 4g scatter (r = −0.581, n = 80) | not reproducible here: the per-organoid NfL values are not in the repository and the regression was run in GraphPad Prism |

### Runtimes

Measured on an 8-core workstation (WSL2) with Python 3.13.13. Every script is
single-threaded (it pins BLAS to one thread). The runtimes in the tables above
are wall-clock times with the machine otherwise idle.
`exp3_perm_power.py` alone took 31.6 min (1,898 s). `run_section.sh --skip-slow` took 10 min (13 steps run).
The complete default run (17 steps, the four slow steps concurrent) was timed twice, both times while
other jobs were occupying the same machine: 63 min and 92 min wall-clock (roughly 50 to 58 min of CPU for
each slow step, against 31 min when one runs alone). Run with `--sequential` on an idle machine it amounts
to about 2 h 20 min: four times 32 min plus the 10 min of the other steps.

Source data for this section is in [`../data/04_drug_response_prediction`](../data/04_drug_response_prediction); see
[`../data/README.md`](../data/README.md) for formats and column schemas.
