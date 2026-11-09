# Chirality peak-ratio analyses re-run with true log10 ratios (2026-11-04)

Change: pairwise features log10(Hi/Hj) with no floor, instead of log10((Hi + 1e-6)/(Hj + 1e-6)).
Fitted peak intensities H are ~1e-14 to 1e-12, so the 1e-6 floor had turned the features into scaled
height differences. All H > 0 in every cache, so no floor is needed. Nothing else changed: same caches,
labels, seeds and pipelines. The original scripts and outputs are untouched.

Scripts here are copies of the originals with only the floor and the output folder changed.
Exception: exp3_peakfit_ratios_truelog.py reads the cached exp3 peak fits (refitting from raw spectra
fails in the current SciPy). With FLO=1e-6 it reproduces the original 0.613 / 0.938 / 0.920 exactly
(log_exp3_peakfit_ratios_checkfloor.txt).

## exp2 (blast classification), exp2_chirality_classify_truelog.py

| | old (+1e-6) | true log ratio |
|---|---|---|
| LOO AUC 66 ratios, 2 h / 24 h | 1.00 / 1.00 | 1.00 / 1.00 |
| Permutation p (B = 5000), z 2 h / 24 h | 2e-4, 5.14 / 5.12 | 2e-4, 5.18 / 5.17 |
| Ratios with q < 0.05, 2 h / 24 h | 13 / 3 of 66 | 45 / 51 of 66 |
| 2 h dilution scores, neg ctrl / 25 / 50 / 75 % | -0.21 / 0.21 / 0.12 / 0.80 | 0.05 / 0.17 / 0.26 / 0.83 |
| 2 h rho_mix (rho_all) | 0.656 (0.869) | 0.781 (0.917) |
| 24 h dilution scores, neg ctrl / 25 / 50 / 75 % | -0.16 / 0.22 / 0.27 / 0.68 | 0.13 / 0.19 / 0.53 / 0.78 |
| 24 h rho_mix (rho_all) | 0.597 (0.887) | 0.870 (0.916) |

## exp3 (drug experiment, 3 h)

| | old (+1e-6) | true log ratio |
|---|---|---|
| LOO AUC 66 ratios (perm_power, PCA15) | 0.938, z = 4.5 | 0.887, z = 4.0 (0/20000 exceed, p < 5e-5) |
| LOO AUC full shape (unchanged) | 0.903, z = 4.3 | 0.903, z = 4.3 |
| 12 area fractions / 78 combined (peakfit_ratios) | 0.613 / 0.920 | 0.613 / 0.895 |
| Ratios with q < 0.05 | 0 of 66 | 17 of 66 |
| Top ratios (Cohen's d) | (8,3)/(10,2) -0.86, (8,3)/(8,4) -0.76, (6,5)/(7,6) +0.77 | (8,3)/(10,2) -1.78, (9,4)/(7,6) +1.39, (10,3)/(10,5) -1.09 |
| Blast vs blast + resveratrol on axis | p = 0.0001 | p = 0.0011 |
| Resveratrol rescue along axis | 43 % (95 % CI 28-60) | 31 % (15-46) |
| MG132 move along axis | 5 % (-10 to 18) | -3 % (-22 to 13) |
| DiD on axis: resveratrol / MG132 | p = 0.48 / 0.0008 | p = 0.79 / 0.023 |
| Organoid score LOO AUC (sham vs blast) | 0.949 | 0.913 |
| Score fold-change: blast+Res / blast+MG132 / sham+Res / sham+MG132 | 1.55 / 1.96 / 0.60 / 1.26 | 1.70 / 2.08 / 0.64 / 1.33 |
| WB correlation, NfL: Pearson r (p), Spearman rho (p) | -0.81 (0.051), -0.83 (0.042) | -0.75 (0.084), -0.77 (0.072) |
| WB correlation, NfH | -0.51 (0.30), -0.54 (0.27) | -0.44 (0.39), -0.37 (0.47) |

## Files
- exp2_chirality_truelog_results.txt, exp2_chirality_truelog.{png,svg}
- log_exp3_drug_axis.txt, exp3_3H_drug_axis.png, svg/exp3_3H_drug_axis__*.svg
- log_exp3_perm_power.txt, exp3_3H_perm_power.png, svg/exp3_3H_perm_power__*.svg
- log_exp3_3H_organoid_score.txt, exp3_3H_organoid_scores.csv
- svg/wb_swcnt_corr__*.svg (plot_wb_swcnt_corr_truelog.py)
- log_exp3_peakfit_ratios.txt, exp3_3H_peakfit_ratios.png
- run_all.sh: how everything was run
