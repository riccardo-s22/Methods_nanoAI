# Data

Source data for the analyses in this repository, in subdirectories matching the
four Results sections. Each section's code README says which of these files its
scripts read.

| Directory | Contents | Size |
|---|---|---|
| `01_swcnt_blast_classification/` | NIR excitation–emission spectra, blast-classification experiment | 5.9 MB |
| `02_geneformer_in_silico_screen/` | in-silico perturbation statistics, pathway and drug-signature enrichment tables | 24 MB |
| `03_protein_corona_proteomics/` | Spectronaut DIA report | 1.7 MB |
| `04_drug_response_prediction/` | NIR spectra of the drug experiment, chirality peak features, per-well scores | 13 MB |

---

## 01_swcnt_blast_classification/

### spectra/

`2h/` and `24h/` are the two reads of the **same** 96-well plate (the 2 h files
are named `Experiment_2_TBI_single__*`, the 24 h files `Experiment_2_TBI_24H__*`).
One file per excitation wavelength: 570, 640, 670, 750 and 780 nm.

Format: 512 emission rows × 97 columns. Column 0 is the emission axis in nm;
columns 1–96 are the wells, row-major, so well index *i* (0-based) is plate row
`i // 12` (A–H) and column `i % 12 + 1`. **Well A1 (index 0) is a synthetic
0→1000 ramp placeholder and is excluded by the loader**, not real data.

Plate layout: sham = columns 2, 4, 6, 8; blast = columns 3, 5, 7, 9 —
interleaved, so a smooth plate gradient is balanced between the groups by
design. Columns 10–12 hold the dilution series and the
SWCNT-in-unconditioned-medium negative control. With `MIX_INVERTED=1`, column 10
is 25% blast and column 12 is 75% blast.

n = 32 sham + 32 blast per read. One organoid = one well = one experimental unit.

### chirality_features/

Cached split-Gaussian fits of the 12 canonical (n,m) species and the permutation
nulls behind Supplementary Figure 3:

- `exp2_{2h,24h}_mix_peakfeats.npz` — per-well peak features
- `exp2_{2h,24h}_permnull.npz` — leave-one-out LDA AUC null, 5,000 label permutations

---

## 02_geneformer_in_silico_screen/

### in_silico_perturbation/ft_excitatory/

Per-gene statistics for the four FT arms. `rescue_*` perturb TBI nuclei toward
the control state; `induction_*` perturb control nuclei toward the TBI state.

| Arm | Genes tested | In the coverage-gated table |
|---|---|---|
| `rescue_delete` | 17,189 | 12,689 |
| `rescue_overexpress` | 17,189 | 12,689 |
| `induction_delete` | 17,223 | 13,393 |
| `induction_overexpress` | 17,222 | 13,392 |

Two tables per arm:

- `*_unfiltered_all_gene_ranked.tsv` — every gene tested. This is the complete
  in-silico perturbation statistics.
- `*_filtered_high_confidence_gene_ranked.tsv` — the subset passing the
  **coverage gate** (perturbed in ≥25 nuclei **or** ≥2% of start-state nuclei),
  minus universal abundance artifacts (MT-\*, RPL\*/RPS\*, MALAT1, hemoglobins).

> **This coverage gate is not the eligibility rule described in the Methods.**
> The paper's gene eligibility is the Geneformer-preprint rule — detected in
> ≥20% of either group's excitatory nuclei, giving the 7,374-gene universe — and
> it was applied **downstream**, in the pathway and drug analyses, not in these
> tables. The absence of `det_start`/`det_goal` columns here is how you can tell
> the two apart. Do not read `filtered_high_confidence` as the paper's eligible
> set, and do not read its row count as the 7,374-gene universe.

Columns:

| Column | Meaning |
|---|---|
| `rank` | rank by mean shift within the arm |
| `gene_symbol`, `ensembl_id`, `gene_token` | gene identifiers |
| `n_cells_perturbed`, `percent_TBI_cells_perturbed` | coverage |
| `mean_shift` | effect size — mean shift toward the goal state (`Shift_to_goal_end`) |
| `median_shift`, `shift_sd`, `shift_sem` | shift distribution |
| `n_donors_contributing` | donors with ≥1 perturbed nucleus |
| `donor_consistency_score` | fraction of donors whose donor-mean shift matches the overall sign |
| `goal_end_vs_random_pval`, `goal_end_fdr`, `sig` | two-sided Wilcoxon rank-sum against the internal random-gene background, Benjamini–Hochberg corrected |

`donor_consistency_score` is meaningful here because the FT split held out three
donors. It is **not** available for the MT screen, which held out one donor per
condition — see below.

**The MT model's in-silico perturbation tables are not included because they no
longer exist.** They lived under `Geneformer_JointPG/results/isp/`, in one of the
two cluster project directories that were deleted; no copy survives. The MT
model's *derived* outputs do survive and are included — its pathway-enrichment
matrix and its drug-candidate tables below.

### pathway_enrichment/

Source matrices for Supplementary Figure 6: `excit_pathway_enrichment_matrix_FT.tsv`
and `..._MT.tsv`, holding −log10 FDR-corrected p per term per arm. Enrichr
(GSEApy), one-sided Fisher exact, FDR-corrected within each library, against GO
Biological Process 2026, KEGG 2026, Reactome 2024 and MSigDB Hallmark 2020, with
the 7,374-gene excitatory background.

The size difference between the two files is real, not truncation: the FT model
tested 7,837 terms with 367 reaching FDR-corrected p < 0.05 in at least one arm,
354 of them in the broad induction/overexpression arm alone.

### drug_signature_enrichment/

`FT/` and `MT/` each hold `drug_candidates_ranked*.tsv` (`_union` = the full
target sets, `_core` = the bidirectional core) and `all_sig_drug_terms.tsv`.
`crossmodel_consensus.tsv` is the intersection, and the source of **Figure 2b**.

Two statistics in that file are easy to confuse: `combined_signal` is the sum of
−log10 adjusted p across signatures and is what compounds were **ranked** by;
`best_adjp_FT` / `best_adjp_MT` is the best single signature per model and is
what Figure 2b **plots**. `in_FT_core` marks membership of the FT bidirectional
core (bold in the figure).

> `crossmodel_consensus.tsv` has **eight** rows, but the paper reports **six**
> consensus drugs. The two extra rows, `vx` and `sodium`, are fragments of
> compound names that the library's name parser split — not compounds. They were
> discarded. The six real consensus drugs are resveratrol, diclofenac,
> vemurafenib, dexamethasone, morphine and methamphetamine.

---

## 03_protein_corona_proteomics/

`Riccado_Report.csv` — the Spectronaut 17 directDIA protein-group report, the
input to `../../03_protein_corona_proteomics/ms_analysis.py`. 6,198 protein
groups; 5,850 quantified in all six analysed samples.

Quantities are in `<run>.raw.PG.Quantity` columns. The report contains more runs
than the paper analyses, including `R`-suffixed re-runs. The six analysed
samples are:

| Group | Samples |
|---|---|
| Control | C4R, C5R, C9 |
| Blast | A1, A2, C6R |

**C4 and C5 are excluded** — they belong to a separate batch; C4R and C5R are
their re-runs. Any other run column in the file is not part of the published
comparison.

**The raw DIA files are not in this repository.** There are 18 of them totalling
3.4 GB, each 167–212 MB, so every one exceeds GitHub's 100 MB per-file limit.
Mass-spectrometry raw data belongs in a ProteomeXchange repository (PRIDE) in any
case, and that is where it should be deposited for the paper; this directory
holds the processed report only.

---

## 04_drug_response_prediction/

### spectra/

Two plates — one sham, one blast — each read at two time points after SWCNT
addition, at five excitation wavelengths: 20 files.

| Files | Read |
|---|---|
| `SHAM_EX_3_3H__*`, `TBI_EX3_3H__*` | 3 h — **the read analysed in the paper** |
| `SHAM_EX3_6H__*`, `TBI_EX3_6H__*` | 6 h — a second read of the same physical wells, not reported |

Format: 512 emission rows × 98 columns. Column 0 is the emission axis; column 1
is a synthetic ramp marker and is excluded; the remaining columns are the wells.

Plate layout, as encoded by `exp3_load.py`:

| Plate | Columns 1–4 | Columns 5–8 | Columns 9–12 |
|---|---|---|---|
| Sham | sham vehicle | resveratrol | MG132 |
| Blast | MG132 | blast vehicle | resveratrol |

32 wells per condition except MG132-on-sham, which is 24 (rows A–F only). Wells
G9–12 and H9–10 are SWCNT in unconditioned medium; H11–12 are medium without
SWCNTs. 184 wells in total.

**The blast plate's resveratrol and MG132 columns are swapped relative to the
plate README shipped with the raw files** (`../../04_drug_response_prediction/README_plate_layout_original.txt`).
The layout above is the corrected one, confirmed against lab records. Use
`exp3_load.py`, not the original README, to assign conditions.

### chirality_features/

- `exp3_3H_peakfeats.npz` — **the authoritative cached split-Gaussian fits.** The
  pipeline reads this rather than refitting: refitting from the raw spectra fails
  under the installed SciPy, so this cache, not the `.txt` files, is the
  effective feature source for this section.
- `exp3_{3H,6H}_chirality_peak_features.csv` — the same features as a flat table,
  one row per well: group, plate row and column, 12 `gauss_max`, 12 `gauss_auc`,
  12 area fractions and the 66 pairwise log10 intensity ratios.

`gauss_auc` falls back to trapezoidal integration when a Gaussian fit fails, so a
few area values are integrals rather than fitted areas. The 66 pairwise log
ratios of 12 quantities have rank 11.

### per_well_scores/

Per-organoid disease-axis scores, one row per well. Sham and blast wells are
scored leave-one-out (standardization and LDA refit without the scored well);
drug-treated wells use the model fitted on all sham and blast wells. Scores are
rescaled so sham = 1 and blast = 2.

| File | Feature definition |
|---|---|
| `exp3_3H_organoid_scores.csv` | `log10((H_i + 1e-6)/(H_j + 1e-6))` — **the published values** |
| `exp3_3H_organoid_scores_truelog.csv` | `log10(H_i/H_j)`, no floor — the sensitivity re-run |
| `exp3_3H_organoid_scores_shape.csv` | full 2,560-dimensional spectral shape instead of peak features |

See the repository README on the floor: because fitted peak intensities are of
order 1e-14 to 1e-12, the 1e-6 floor dominates the ratio and turns these features
into scaled height differences. The main text reports the floored values; the
other two files let the dependency be inspected.

---

## Not included

- **DIA raw files** — 3.4 GB, over GitHub's per-file limit; belongs in PRIDE.
  See section 03 above.
- **snRNA-seq data** — public, at NCBI GEO
  [GSE209552](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE209552)
  (Garza et al. 2023). The preparation scripts in
  `../02_geneformer_in_silico_screen/01_dataset_preparation/` download and
  process it.
- **Model weights** — Geneformer-V2-104M is at
  [ctheodoris/Geneformer](https://huggingface.co/ctheodoris/Geneformer). The
  fine-tuned FT checkpoint is available from the corresponding author; the MT
  checkpoint was deleted and no copy survives.
- **MT in-silico perturbation tables** — deleted, see section 02 above.
- **Western blot images and densitometry** — quantified in ImageJ, plotted in
  GraphPad Prism; available from the corresponding author.
- **7-day spectral read** — a separate experiment, not reported in the paper.
