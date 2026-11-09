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

Environment: `cwd = exp2/`, with `EXP2_DIR=exp2/{2h,24h}`, `EXP2_TAG={2h,24h}`,
`MIX_INVERTED=1`.

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

## Inputs

Five `.txt` files per time point, one per excitation wavelength (570, 640, 670,
750, 780 nm), each 512 emission rows × 97 columns (column 0 = emission axis,
columns 1–96 = the wells row-major). Well A1 is a synthetic 0→1000 ramp
placeholder and is excluded by the loader.
