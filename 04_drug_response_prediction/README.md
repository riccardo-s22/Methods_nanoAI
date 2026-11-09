# Results 4 — early prediction of a drug's neuroprotective effect

**Figure 4; Supplementary Figures 7–10.**

Organoids were blast- or sham-exposed, treated for 12 h with 100 µM resveratrol
(nominated by the Geneformer screen), 10 µM MG132 or vehicle, and their
conditioned medium read with the SWCNT sensor 36 h after injury. The same
organoids were then collected on day 7 for western blotting, so the optical
readout can be tested against the biochemical outcome it is meant to anticipate.

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

Environment: `cwd` = the directory holding `exp3_3H_peakfeats.npz` (the cache
with corrected column labels), with `EXP3_TP=3H`.

## figures/

`plot_wb_swcnt_corr.py` — **Supplementary Figure 10b–c**. Correlates the mean SWCNT score of each of the six conditions (sham or blast ×
vehicle, resveratrol, MG132) against the mean NfL and NfH western-blot values of
the same condition, normalized to GAPDH and expressed relative to vehicle-treated
sham, by Pearson and Spearman (n = 6 conditions). Data: `corr_wb_swcnt.xlsx`.

Figure 4a (timeline), 4c (representative blots) and 4d–f (day-7 p-tau, NfH and
NfL quantification), together with Supplementary Figures 9 and 10a, involve no
code: blots were quantified in ImageJ and the bar charts drawn in GraphPad Prism
(two-way ANOVA with Tukey's multiple-comparisons test).

## supporting_checks/

| Script | Question |
|---|---|
| `exp3_3H_pressure.py` | Does the headline sham-versus-blast result survive the spatial plate gradient and outlier wells? Unlike the section 1 plate, sham and blast here sit on different plates, so position is not balanced by design and has to be checked. |
| `exp3_drug_swap.py` | Tests both candidate labellings of the blast plate's drug columns. The sham→blast axis is unchanged by the swap, since it is built from sham and blast wells, neither of which is affected. |

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

## Inputs

| File | What it is |
|---|---|
| `SHAM_EX_3_3H__{570,640,670,750,780}.txt`, blast equivalents | Raw NIR emission, one file per excitation per plate per time point. |
| `exp3_3H_peakfeats.npz` | Cached split-Gaussian fits with corrected labels. |
| `corr_wb_swcnt.xlsx` | Western-blot densitometry for the SWCNT-versus-biochemistry correlation. |

## Known limitations

- **The peak fitting is not currently reproducible from raw spectra.** The
  pipeline loads the cached `exp3_3H_peakfeats.npz`; refitting the split
  Gaussians from the raw `.txt` files fails under the installed SciPy. The cache
  is the authoritative feature matrix for this section.
- `gauss_auc` falls back to trapezoidal integration when a Gaussian fit fails,
  so a small number of area values are not fitted areas.
- The 66 pairwise log ratios of 12 quantities have rank 11, so a PCA asking for
  more than 11 components returns noise directions.

Source data for this section is in [`../data/04_drug_response_prediction`](../data/04_drug_response_prediction); see
[`../data/README.md`](../data/README.md) for formats and column schemas.
