# Results 1 — SWCNT nano-sensor classifies the blast TBI cortical organoid model

**Figure 1; Supplementary Figures 1–4.**

(GT)15-SWCNT sensors were incubated for 24 h in medium conditioned by
blast-exposed or sham cortical organoids, read at five excitation wavelengths,
and classified. Trained classifiers were then applied without retraining to
three dilutions of blast- in sham-conditioned medium (25%, 50%, 75%).

## Design

One organoid = one well = one experimental unit. n = 32 sham + 32 blast per time
point, on a single 96-well plate with **interleaved** columns (sham = columns
2, 4, 6, 8; blast = 3, 5, 7, 9), so any smooth plate gradient is balanced
between the groups by design rather than by correction. Columns 10–12 carry the
dilution series and the SWCNT-in-unconditioned-medium control.

Two reads of the same plate, 2 h and 24 h. The 24 h read is Figure 1; the 2 h
read is Supplementary Figure 1. Each well is a 2,560-dimensional shape vector
(5 excitations × 512 emission channels, L1-normalized per excitation to remove
well-to-well brightness differences).

## Run order

| Script | What it does |
|---|---|
| `exp2_load.py` | Loader. Parses the five per-excitation `.txt` files into the well × 2,560 matrix and assigns the plate layout. Imported by everything else; run `MIX_INVERTED=1` so that column 10 = 25% blast and column 12 = 75% blast. |
| `exp2_analysis.py` | Leave-one-out classification on raw and shape-normalized spectra; PCA; PERMANOVA of sham-versus-blast separation. |
| `exp2_analysis2.py` | **Figure 1a** (mean spectra at ex 670 nm), **Figure 1b** (shape-spectrum PCA) and **Supplementary Figure 1**. Also the permutation null for the LDA, per-excitation separability, and the broad-dimming and peak-shift decomposition. |
| `exp2_classifier_compare.py` | **Figure 1d**. Dilution dose-response across the five classifiers (LDA, L2 logistic regression, linear SVM, elastic-net logistic regression on 20 PCs, elastic-net on all 2,560 features). Each score is rescaled to sham mean 0, blast mean 1 for a common dose axis. |
| `exp2_dose_pca_figs.py` | LDA decision-score dose-response and the PCA of all wells showing the dilutions nesting between sham and blast. |
| `exp2_mixes.py` | **Figure 1c**, with `exp2_mix_intermediate.py`. Projects the dilution series and the negative control onto the sham→blast axis and builds the difference-from-sham spectra of blast medium, the 25/50/75% dilutions and the negative control. |
| `exp2_mix_intermediate.py` | Contributes the per-group difference-from-sham spectra of **Figure 1c**. Tests whether the dilutions are a genuinely intermediate phenotype without relying on the whitened LDA axis: un-whitened centroid projection, cosine similarity of each group's difference spectrum to the blast difference spectrum, and projection magnitude. |
| `exp2_mix_features.py` | Per-feature discrimination AUC against dose-monotonicity, to find features that are both discriminative and dose-graded. |
| `exp2_2h_vs_24h.py` | Whether the 24 h signature is the matured 2 h signature or a different spectral change. |

Environment: no working directory is assumed. `EXP2_DIR` selects the spectra folder of the
read, `EXP2_TAG={2h,24h}` labels titles and output file names (set the two together), and
`MIX_INVERTED=1` gives the published dilution labelling. The exact commands are in
[Running this section](#running-this-section).

## figures/

| Script | Output |
|---|---|
| `make_supp_auc_compute.py` → `make_supp_auc_figure.py` | **Supplementary Figure 2**. ROC curves, held-out scores, permutation null, AUC by held-out unit. Two steps: compute writes `supp_auc_results.json`, figure draws it. |
| `exp2_chirality_classify.py` | **Supplementary Figure 3**. Chirality-resolved peak features — split-Gaussian fits of the 12 canonical (n,m) species per well — and the permutation test of sham-versus-blast separation on them. |
| `make_suppfig4_dilution_readouts.py`, `make_suppfig4_gaussian_readouts.py` | **Supplementary Figure 4**. The two readouts of the 24 h dilution series. |
| `exp2_pca_quant.py` | PCA variance quantification reported in the text. |

`make_supp_auc_compute.py` deliberately differs from the original `exp2`
pipeline in three ways, all documented in
`supporting_checks/README_supp_auc_original.txt`: `PCA(20, random_state=0)`
instead of `PCA(20)` (on a 64 × 2560 matrix scikit-learn silently selects the
randomized SVD solver, so an unseeded subspace differs between runs); AUC from
`decision_function` rather than `predict_proba` (LDA probabilities saturate at
exactly 1.0, so wells tie at the float ceiling and `roc_auc_score` scores those
pairs 0.5); and scaler plus PCA refitted **inside** every CV fold. The first two
together made the 2 h held-out AUC print as either 1.000 or 0.9844 depending on
the run.

## supporting_checks/

Not figures. These underwrite claims the paper makes about the method.

| Script | Question |
|---|---|
| `check_simple_descriptors.py` | Could simple univariate descriptors replace the classifier? At 24 h largely yes — the emission centroid red-shifts 6.0–8.5 nm at every excitation, best single-descriptor AUC 0.858. At 2 h no — best centroid AUC 0.659 (p = 0.145 after correcting for five excitations), total intensity AUC 0.588. The 2 h signal is multivariate spectral shape with no univariate correlate, so the classifier is load-bearing there, not decorative. |
| `check_permutation_power.py` | Power of the permutation tests, and how many permutations are needed to certify significance. |

## sensitivity_true_log_ratios/

`exp2_chirality_classify_truelog.py` re-runs Supplementary Figure 3 with true
`log10(H_i/H_j)` ratios instead of the floored form — see the repository README.
The classification result is unchanged (LOO AUC 1.00 at both time points); what
changes is the number of individually significant ratios (13→45 of 66 at 2 h,
3→51 at 24 h) and the dilution-series correlations, which improve.

## Running this section

Everything runs on an ordinary CPU machine from the data shipped in
[`../data/01_swcnt_blast_classification`](../data/01_swcnt_blast_classification), from any working
directory, with no path to edit:

    bash 01_swcnt_blast_classification/run_section.sh               # everything, about 65 minutes
    bash 01_swcnt_blast_classification/run_section.sh --skip-slow   # the quick steps only, about 2 minutes
    bash 01_swcnt_blast_classification/run_section.sh --skip-svg    # without the SVG re-rendering
    bash 01_swcnt_blast_classification/run_section.sh --only '_24h$' # only the steps whose name matches a regular expression

The runner executes the scripts in the order of the tables above, once per read (2 h and 24 h)
where a script analyses one read, writes every output under
`results/01_swcnt_blast_classification/`, keeps each script's stdout and stderr in
`results/01_swcnt_blast_classification/logs/<step>.log`, prints one `PASS` / `FAIL` / `SKIP` line
per step and exits non-zero if any step failed. `results/` is not tracked by git; `data/` is never
written to.

Requirements: Python 3 with NumPy, SciPy, scikit-learn, pandas, matplotlib, joblib and Pillow
(versions used: Python 3.13, NumPy 2.5.1, SciPy 1.18.0, scikit-learn 1.9.0, pandas 3.0.3,
matplotlib 3.11.0, joblib 1.5.3). No GPU, no network access.

### Environment variables

| Variable | Meaning | Default |
|---|---|---|
| `EXP2_TAG` | Label of the read, `2h` or `24h`, used in figure titles and output file names. It does **not** select the spectra: that is `EXP2_DIR` | `2h` |
| `EXP2_DIR` | Folder holding the five spectra `.txt` files of one read. This is what selects the read | `data/01_swcnt_blast_classification/spectra/2h` |
| `MIX_INVERTED` | Set (to `1`) for the published labelling, column 10 = 25% blast and column 12 = 75% blast. Every script that uses the dilution wells sets it itself if it is unset | unset in `exp2_load.py` |
| `EXP2_OUT` | Output folder of the top-level scripts | `results/01_swcnt_blast_classification` |
| `EXP2_PEAKFEATS_DIR` | Folder holding `exp2_{2h,24h}_peakfeats.npz` (see *Input not in the repository*) | `data/01_swcnt_blast_classification/chirality_features` |
| `FEAT` | `make_suppfig4_gaussian_readouts.py` only: `paper` (floored ratios) or `truelog` | `paper` |
| `NANOAI_DATA`, `NANOAI_RESULTS` | Roots of the input and output trees, for every script of the repository | `<repo>/data`, `<repo>/results` |

With no variable set, each script runs its original default: the top-level scripts analyse the
**2 h** read; the scripts in `figures/`, `supporting_checks/` and `sensitivity_true_log_ratios/`
analyse both reads, or the 24 h read where they are written for it. The published Figure 1 is the
24 h read, so the top-level scripts need **both** `EXP2_TAG=24h` and
`EXP2_DIR=data/01_swcnt_blast_classification/spectra/24h` for it (the runner does both reads).
As in the original scripts, `EXP2_TAG` on its own only changes the label: `EXP2_TAG=24h` without
`EXP2_DIR` analyses the 2 h spectra and names the outputs `24h`.

### Commands, inputs, outputs and run times

Commands are given from the repository root; `$TP` is `2h` or `24h` (set it first, for example
`TP=24h`), and `$SP` stands for `data/01_swcnt_blast_classification/spectra` (`SP=data/01_swcnt_blast_classification/spectra`).
Inputs are relative to
`data/01_swcnt_blast_classification/`, outputs to `results/01_swcnt_blast_classification/`.
Run times are wall-clock seconds measured with the runner on the machine used for the paper
(8-core workstation, WSL2), 2 h read / 24 h read where a script runs once per read. The steps
that take more than a minute were timed with the machine otherwise idle, except
`make_supp_auc_compute.py`, which shared the machine with other jobs for part of its run; on a
busy machine the same steps took two to five times longer. Of the 65 minutes of a full run, 40
are `make_supp_auc_compute.py` and about 10 the SVG re-rendering.

| Script | Command | Reads | Writes | Run time |
|---|---|---|---|---|
| `exp2_load.py` | `EXP2_TAG=$TP EXP2_DIR=$SP/$TP MIX_INVERTED=1 python3 01_swcnt_blast_classification/exp2_load.py` | `spectra/$TP/` | nothing (prints the plate layout check) | 1 / 1 s |
| `exp2_analysis.py` | `EXP2_TAG=$TP EXP2_DIR=$SP/$TP OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 01_swcnt_blast_classification/exp2_analysis.py` | `spectra/$TP/` | nothing (statistics on stdout) | 68 / 76 s |
| `exp2_analysis2.py` | `EXP2_TAG=$TP EXP2_DIR=$SP/$TP python3 01_swcnt_blast_classification/exp2_analysis2.py` | `spectra/$TP/` | `exp2_$TP_sham_tbi.png` | 86 / 102 s |
| `exp2_classifier_compare.py` | `EXP2_TAG=$TP EXP2_DIR=$SP/$TP MIX_INVERTED=1 python3 01_swcnt_blast_classification/exp2_classifier_compare.py` | `spectra/$TP/` | `exp2_$TP_classifier_dose_compare.png` | 192 / 226 s |
| `exp2_dose_pca_figs.py` | `EXP2_TAG=$TP EXP2_DIR=$SP/$TP MIX_INVERTED=1 python3 01_swcnt_blast_classification/exp2_dose_pca_figs.py` | `spectra/$TP/` | `exp2_$TP_decision_dose_response.png`, `exp2_$TP_pca_dose.png` | 4 / 5 s |
| `exp2_mixes.py` | `EXP2_TAG=$TP EXP2_DIR=$SP/$TP MIX_INVERTED=1 python3 01_swcnt_blast_classification/exp2_mixes.py` | `spectra/$TP/` | `exp2_$TP_dose_drivers.png` | 5 / 3 s |
| `exp2_mix_intermediate.py` | `EXP2_TAG=$TP EXP2_DIR=$SP/$TP MIX_INVERTED=1 python3 01_swcnt_blast_classification/exp2_mix_intermediate.py` | `spectra/$TP/` | `exp2_$TP_mix_intermediate.png` | 3 / 4 s |
| `exp2_mix_features.py` | `EXP2_TAG=$TP EXP2_DIR=$SP/$TP MIX_INVERTED=1 python3 01_swcnt_blast_classification/exp2_mix_features.py` | `spectra/$TP/` | `exp2_$TP_mix_feature_scan.png` | 7 / 4 s |
| `exp2_2h_vs_24h.py` | `python3 01_swcnt_blast_classification/exp2_2h_vs_24h.py` | `spectra/2h/`, `spectra/24h/` | `exp2_2h_vs_24h_diffspectrum.png` | 2 s |
| `figures/make_supp_auc_compute.py` | `python3 01_swcnt_blast_classification/figures/make_supp_auc_compute.py` | `spectra/2h/`, `spectra/24h/` | `supp_fig2_auc/supp_auc_results.json` | 2400 s |
| `figures/make_supp_auc_figure.py` | `python3 01_swcnt_blast_classification/figures/make_supp_auc_figure.py` | `results/…/supp_fig2_auc/supp_auc_results.json` (the step above) | `supp_fig2_auc/supp_auc_blast_vs_sham.{svg,pdf,png}` | 4 s |
| `figures/exp2_chirality_classify.py` | `python3 01_swcnt_blast_classification/figures/exp2_chirality_classify.py` | `chirality_features/exp2_{2h,24h}_peakfeats.npz` (**not shipped**, see below), `chirality_features/exp2_{2h,24h}_mix_peakfeats.npz`, `chirality_features/exp2_{2h,24h}_permnull.npz` | `supp_fig3_chirality/exp2_chirality_results.txt`, `supp_fig3_chirality/exp2_chirality.{png,svg}` | 6 s |
| `figures/make_suppfig4_dilution_readouts.py` | `python3 01_swcnt_blast_classification/figures/make_suppfig4_dilution_readouts.py` | `spectra/24h/` | `supp_fig4_dilution/suppfig4_dilution_readouts.{svg,png,tif}` | 4 s |
| `figures/make_suppfig4_gaussian_readouts.py` | `python3 01_swcnt_blast_classification/figures/make_suppfig4_gaussian_readouts.py` | `chirality_features/exp2_24h_peakfeats.npz` (**not shipped**), `chirality_features/exp2_24h_mix_peakfeats.npz` | `supp_fig4_dilution/suppfig4_gaussian_readouts.{svg,png,tif}` | 4 s |
| `figures/exp2_pca_quant.py` | `python3 01_swcnt_blast_classification/figures/exp2_pca_quant.py` | `spectra/2h/`, `spectra/24h/` | `pca_quant/exp2_pca_quant_results.txt` | 13 s |
| `supporting_checks/check_simple_descriptors.py` | `python3 01_swcnt_blast_classification/supporting_checks/check_simple_descriptors.py` | `spectra/2h/`, `spectra/24h/` | nothing (statistics on stdout) | 2 s |
| `supporting_checks/check_permutation_power.py` | `python3 01_swcnt_blast_classification/supporting_checks/check_permutation_power.py` | `spectra/2h/`, `spectra/24h/` | nothing (statistics on stdout) | 3 s |
| `sensitivity_true_log_ratios/exp2_chirality_classify_truelog.py` | `python3 01_swcnt_blast_classification/sensitivity_true_log_ratios/exp2_chirality_classify_truelog.py` | `chirality_features/exp2_{2h,24h}_peakfeats.npz` (**not shipped**), `chirality_features/exp2_{2h,24h}_mix_peakfeats.npz` | `sensitivity_true_log_ratios/exp2_chirality_truelog_results.txt`, `sensitivity_true_log_ratios/exp2_chirality_truelog.{png,svg}`, `sensitivity_true_log_ratios/exp2_{2h,24h}_permnull.npz` | 369 s on the first run, seconds afterwards |
| `figures/make_suppfig4_gaussian_readouts.py`, true log ratios | `FEAT=truelog python3 01_swcnt_blast_classification/figures/make_suppfig4_gaussian_readouts.py` | as above | `supp_fig4_dilution/suppfig4_gaussian_readouts_truelog.{svg,png,tif}` | 2 s |

Scripts that print their statistics only (`exp2_analysis.py`, the two checks) leave them in the
runner's log, `results/01_swcnt_blast_classification/logs/<step>.log`.

`exp2_analysis.py` is the only script of this section that does not pin the BLAS thread count
itself. Run it with one thread, as above and as the runner does: on these 64 × 2,560 matrices
an eight-thread BLAS was far slower here (stopped unfinished after more than 30 minutes and 60
CPU-minutes on a busy machine, against about 2 CPU-minutes with one thread), and the printed
statistics with one thread are identical to the original logs.

`exp2_chirality_classify.py` reads the shipped permutation null instead of recomputing it. The
true-log-ratio script does not: the shipped null belongs to the floored features, so it computes
its own 5,000-permutation null on the first run (eight joblib workers) and caches it with its
outputs. Both scripts read the shipped dilution-well fits; if those were missing they would refit
them with `vendored/chirality_interpolated_extraction.py` into `results/`.

### SVG panels

The panels of Figure 1 and Supplementary Figure 1 were exported as SVG with
[`../shared/make_svg.py`](../shared), which runs a figure script unchanged with an SVG backend and
splits the figure into its panels:

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 shared/make_svg.py 01_swcnt_blast_classification/exp2_analysis2.py 01_swcnt_blast_classification/svg/exp2_24h EXP2_TAG=24h EXP2_DIR=data/01_swcnt_blast_classification/spectra/24h MIX_INVERTED=1

Set the three thread variables in the shell, as shown: the figure scripts pin the BLAS thread
count themselves, but under `make_svg.py` NumPy is already loaded when they do, so the pin comes
too late, and with a multi-threaded BLAS `exp2_analysis2.py` and `exp2_classifier_compare.py` run
many times slower (the command above was stopped unfinished after 10 minutes without them).

The runner does this for `exp2_analysis2.py`, `exp2_classifier_compare.py`,
`exp2_dose_pca_figs.py`, `exp2_mixes.py` and `exp2_mix_intermediate.py` at both reads, into
`results/01_swcnt_blast_classification/svg/exp2_{2h,24h}/` (`*__full.svg`, `*__panel_X.svg`).
Each re-rendering runs the script again, so it takes as long as the script itself:
98 / 96 s for `exp2_analysis2.py`,
213 / 203 s for
`exp2_classifier_compare.py`, seconds for the other three.

### Input not in the repository

`figures/exp2_chirality_classify.py`, `figures/make_suppfig4_gaussian_readouts.py` and
`sensitivity_true_log_ratios/exp2_chirality_classify_truelog.py` read the cached split-Gaussian
fits of the 32 sham and 32 blast wells, `exp2_2h_peakfeats.npz` and `exp2_24h_peakfeats.npz`
(arrays `MXs`, `AUs`, `MXt`, `AUt`, each 32 wells × 12 chiralities). These two files are **not**
in `data/01_swcnt_blast_classification/chirality_features/`, which holds only the dilution-well
fits (`*_mix_peakfeats.npz`) and the permutation nulls (`*_permnull.npz`). Until they are added,
the runner prints `SKIP` for these three scripts (Supplementary Figure 3, the Gaussian version of
Supplementary Figure 4 and the true-log-ratio re-run). To run them, put the two files in that
folder or point `EXP2_PEAKFEATS_DIR` at the folder that holds them.

### vendored/

`vendored/chirality_interpolated_extraction.py` is the split-Gaussian extractor (peak
localization, windowing and fitting), copied unmodified. In this section it is only imported when
the dilution-well fits have to be rebuilt, which does not happen while the shipped caches are
present.

### What a rerun reproduces

Checked on the machine used for the paper, with the package versions listed above:

- The statistics printed by every top-level script are identical, at both reads, to the logs of
  the original runs of these seeded scripts.
- `pca_quant/exp2_pca_quant_results.txt` is identical to the original file.
- `supp_fig2_auc/supp_auc_results.json`: the 30 AUCs, the two 1,000-permutation nulls, the LDA
  held-out scores and all ROC curves are identical to the original JSON. The solver of the two
  elastic-net models and `SVC` are not seeded in `make_supp_auc_compute.py`; between runs the
  elastic-net score gaps (`gap`, values of 3 to 8) moved by up to 0.2 and nothing else changed.
  The figure does not use them: `supp_auc_blast_vs_sham.png` is identical to the original.
- `figures/make_suppfig4_dilution_readouts.py` uses an unseeded `PCA(20)` (randomized solver);
  its printed scores, projections and cosines were the same in repeated runs, but the image is
  not identical from run to run (point positions in panel a move by a fraction of a point).
- With the two `peakfeats` files (see *Input not in the repository*):
  `supp_fig3_chirality/exp2_chirality_results.txt` and
  `sensitivity_true_log_ratios/exp2_chirality_truelog_results.txt` are identical to the original
  files apart from the lines that print paths and the cache-reload notice; the recomputed
  true-log-ratio nulls are identical to the original ones; and with the shipped caches removed,
  the dilution-well fits and the 5,000-permutation nulls are rebuilt identical to the shipped
  `*_mix_peakfeats.npz` and `*_permnull.npz`.
- The SVG panels of `exp2_analysis2.py` and `exp2_mix_intermediate.py`, at both reads, are
  identical to the originally exported SVG files apart from element identifiers, the creation
  date and the last digit of a few coordinates. Those of `exp2_classifier_compare.py` and
  `exp2_dose_pca_figs.py` are not: the original exports were made before the seeds these two
  scripts now carry (solver seed, point jitter), so marker positions and a few plotted values
  differ slightly.

`make_supp_auc_figure.py` asks for the Arial font; where it is not installed Matplotlib falls
back to its default font and prints a warning.

## Inputs

Five `.txt` files per time point, one per excitation wavelength (570, 640, 670,
750, 780 nm), each 512 emission rows × 97 columns (column 0 = emission axis,
columns 1–96 = the wells row-major). Well A1 is a synthetic 0→1000 ramp
placeholder and is excluded by the loader.

Source data for this section is in [`../data/01_swcnt_blast_classification`](../data/01_swcnt_blast_classification); see
[`../data/README.md`](../data/README.md) for formats and column schemas.
