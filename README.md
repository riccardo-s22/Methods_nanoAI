# Coupling foundation AI with nanosensor profiling of organoid culture media for neuroprotective drug discovery

Analysis code for Sirtori, Monroy Lopez, Roxbury and Fallini.

Ryan Institute for Neuroscience and College of Engineering, University of Rhode Island.
Correspondence: riccardo-s@uri.edu

The pipeline links a Geneformer in-silico perturbation screen of a public human
TBI single-nucleus RNA-seq cohort to (GT)15-ssDNA SWCNT nanosensor readouts of
blast-injured human cortical organoids. This repository holds the scripts that
produced the figures and reported statistics, **organized by Results section**.

---

## Layout

| Directory | Results section | Figures |
|---|---|---|
| [`01_swcnt_blast_classification/`](01_swcnt_blast_classification) | SWCNT nano-sensor classifies blast TBI cortical organoid model | Fig. 1; Supp. Figs. 1–4 |
| [`02_geneformer_in_silico_screen/`](02_geneformer_in_silico_screen) | TBI-induced dysregulation of the mechanical pathway … by a foundational transformer model | Fig. 2; Supp. Figs. 5–6 |
| [`03_protein_corona_proteomics/`](03_protein_corona_proteomics) | SWCNT nano-sensor's corona carries brain-specific markers and shows condition-specific differences | Fig. 3 |
| [`04_drug_response_prediction/`](04_drug_response_prediction) | SWCNT nanosensor enables early prediction of drugs' neuroprotective effect | Fig. 4; Supp. Figs. 7–10 |
| [`shared/`](shared) | — | SVG regeneration harness used by sections 1 and 4 |
| [`data/`](data) | all four | source spectra, proteomics report, Geneformer result tables |

Each section directory has its own README mapping every script to the panel or
statistic it produces, in run order.

Within a section:

- top level — the analysis pipeline, in run order
- `figures/` — scripts whose output is a published panel
- `supporting_checks/` — robustness checks behind a methodological choice stated
  in the paper; not themselves figures
- `sensitivity_true_log_ratios/` — see *Chirality feature definition* below

## Figure index

| Panel | Script |
|---|---|
| Fig. 1a, 1b | `01_swcnt_blast_classification/exp2_analysis2.py` |
| Fig. 1c | `01_swcnt_blast_classification/exp2_mixes.py`, `exp2_mix_intermediate.py` |
| Fig. 1d | `01_swcnt_blast_classification/exp2_classifier_compare.py` |
| Fig. 2a | `02_geneformer_in_silico_screen/figures/mt_cls_umap_RECONSTRUCTED.py` (see note below) |
| Fig. 2b | `02_geneformer_in_silico_screen/figures/make_fig2_panelB.py` |
| Fig. 3 | `03_protein_corona_proteomics/ms_analysis.py` |
| Fig. 4a | no script — experimental timeline |
| Fig. 4b | `04_drug_response_prediction/exp3_drug_axis.py` |
| Fig. 4c | no script — representative western blots |
| Fig. 4d–f | no script — GraphPad bar charts of day-7 densitometry |
| Fig. 4g | `04_drug_response_prediction/exp3_3H_organoid_score.py` (scores; scatter assembled outside the script) |
| Supp. Fig. 1 | `01_swcnt_blast_classification/exp2_analysis2.py` (`EXP2_TAG=2h`) |
| Supp. Fig. 2 | `01_swcnt_blast_classification/figures/make_supp_auc_compute.py` → `make_supp_auc_figure.py` |
| Supp. Fig. 3 | `01_swcnt_blast_classification/figures/exp2_chirality_classify.py` |
| Supp. Fig. 4 | `01_swcnt_blast_classification/figures/make_suppfig4_dilution_readouts.py`, `make_suppfig4_gaussian_readouts.py` |
| Supp. Fig. 5a | `02_geneformer_in_silico_screen/figures/umap_all_celltypes_expression.py` (gene-expression UMAP) |
| Supp. Fig. 5b | `02_geneformer_in_silico_screen/figures/umap_excit_ft.py` (FT CLS embeddings) |
| Supp. Fig. 5, metric panels | `02_geneformer_in_silico_screen/figures/make_suppfig5_model_metrics.py` |
| Supp. Fig. 6 | `02_geneformer_in_silico_screen/04_pathway_enrichment/excit_pathway_heatmap_FT.py`, `excit_pathway_heatmap_MT.py` |
| Supp. Fig. 7 | `04_drug_response_prediction/exp3_perm_power.py` |
| Supp. Fig. 8 | `04_drug_response_prediction/exp3_drug_axis.py` |
| Supp. Fig. 9 | no script — GraphPad bar charts from western-blot densitometry |
| Supp. Fig. 10b–c | `04_drug_response_prediction/figures/plot_wb_swcnt_corr.py` |
| Supp. Fig. 10a | no script — GraphPad |

Western blots and their bar charts were quantified in ImageJ and plotted in
GraphPad Prism; no code is involved and none is included.

---

## Note on the MT scripts

The seven `*_RECONSTRUCTED` scripts (the six in
`02_geneformer_in_silico_screen/.../mt_multitask_RECONSTRUCTED/` and
`02_geneformer_in_silico_screen/figures/mt_cls_umap_RECONSTRUCTED.py`) were
rewritten from the Methods after the original MT code was lost. They document
the procedure but were not used to generate the published results, and the MT
checkpoint they need is not available. All other scripts are the original code.

---

## Chirality feature definition

The published chirality analyses compute pairwise features as
`log10((H_i + 1e-6) / (H_j + 1e-6))`. Fitted peak intensities `H` are of order
1e-14 to 1e-12, so that floor dominates the ratio and turns the features into
scaled height *differences* rather than log ratios.

`sensitivity_true_log_ratios/` re-runs the affected analyses with the true ratio
`log10(H_i / H_j)` and no floor. All `H` are strictly positive in every cache, so
no floor is needed. Nothing else changes — same caches, labels, seeds and
pipelines. `04_drug_response_prediction/sensitivity_true_log_ratios/README_truelog.md`
tabulates old against new for every affected statistic. The blast-classification
result is unchanged (LOO AUC 1.00 at both time points); the drug-experiment
chirality AUC moves from 0.938 to 0.887, and the day-7 NfL correlation loses
nominal significance.

The main text reports the floored values. The sensitivity re-run is included so
the dependency is inspectable.

---

## Data

Source data is in [`data/`](data), in subdirectories matching the four Results
sections — see [`data/README.md`](data/README.md) for formats, plate layouts and
column schemas.

| Directory | Contents |
|---|---|
| `data/01_swcnt_blast_classification/` | NIR excitation–emission spectra (2 h and 24 h reads), chirality peak features and permutation nulls |
| `data/02_geneformer_in_silico_screen/` | FT in-silico perturbation statistics (4 arms), pathway-enrichment matrices (FT and MT), drug-signature enrichment and cross-model consensus tables |
| `data/03_protein_corona_proteomics/` | Spectronaut 17 directDIA protein-group report |
| `data/04_drug_response_prediction/` | NIR spectra of the drug experiment (3 h read), chirality peak feature tables, per-well disease-axis scores |

Not in the repository, and why:

- **DIA raw files** — 18 files, 3.4 GB, each 167–212 MB, so each exceeds GitHub's
  100 MB per-file limit. Mass-spectrometry raw data belongs in a ProteomeXchange
  repository (PRIDE) and should be deposited there for the paper. `data/` holds
  the processed Spectronaut report.
- **snRNA-seq data** — public, at NCBI GEO
  [GSE209552](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE209552)
  (Garza et al. 2023). The scripts in `01_dataset_preparation/` fetch and process it.
- **Geneformer-V2-104M** — [ctheodoris/Geneformer](https://huggingface.co/ctheodoris/Geneformer).
- **MT in-silico perturbation tables and the MT checkpoint** — not available.
  The MT model's derived pathway and drug tables are included.
- **Western blots and densitometry** — ImageJ and GraphPad Prism; available from
  the corresponding author.

Scripts carry the absolute paths of the machines they were run on — a Windows
workstation (`/mnt/c/...`, `/mnt/d/...`) for the spectral and proteomic analyses,
and the URI Unity cluster (`/work/pi_cfallini_uri_edu/...`) for the Geneformer
work. They are left as they were run rather than rewritten into portable form,
so paths need adjusting to point at `data/` before re-use.

## Software

Spectral analyses: Python 3.13, NumPy 2.5, SciPy 1.18, scikit-learn 1.9,
pandas 3.0, matplotlib 3.11. Geneformer work: Python 3.10, scanpy 1.11.5,
anndata 0.11.4, Geneformer V2, PyTorch with bf16 on single NVIDIA A100 80GB
GPUs, Slurm array jobs. All statistical tests are two-sided.

## Not included

This repository holds the scripts behind the paper. Left out, and where they
live in the source tree:

- `exp1/` — a superseded pilot of the blast-classification experiment
  (convolutional autoencoder on 16 vs 16 wells, `dim477` latent features). The
  published analysis is the 32-vs-32 interleaved design in section 1.
- `exp2_unsup*.py`, `exp2_zprime.py`, `exp3_zprime.py`, `exp3_umap.py`,
  `exp3_pairdist.py`, `exp3_sham_mg132sham.py`, `exp3_vs_exp2_chirality.py`,
  the ElasticNet variants and the 6 h and 7 d spectral time points —
  exploratory analyses not reported.
- the zero-shot in-silico perturbation pipeline on the pretrained model
  (`scripts/02`–`11`, `perturb_cells_paperQC_nGene1000*.py`,
  `zeroshot_20pct_top100_enrichr.py`) and its astrocyte, microglia and
  oligodendrocyte arms. The paper's screen uses the two **fine-tuned** models.
- `paper/reviewer_concerns/` — analyses written for a reviewer response, not
  for the paper. Two of its scripts are published as Supp. Fig. 2 and appear
  here under `01_swcnt_blast_classification/`.
