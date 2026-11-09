# Results 2 — Geneformer in-silico perturbation screen

**Figure 2; Supplementary Figures 5–6.**

Two models were fine-tuned from Geneformer-V2-104M on the human cortical
snRNA-seq cohort of Garza et al. 2023 (GEO GSE209552; 12 severe-TBI donors,
3 non-neurological controls, 15,172 nuclei after QC):

- **MT** — multi-task, all six major cortical cell types, three heads on the
  final-layer CLS token: disease state (2 classes), cortical region (3), cell
  type (6), task weights 0.40 / 0.20 / 0.40. Corpus depth-matched to 9,428
  nuclei (4,714 control, 4,714 TBI). Retained checkpoint: step 600 of 1,338,
  mean macro F1 0.721.
- **FT** — binary TBI-versus-control classifier on excitatory neurons only
  (6,025 nuclei), accuracy 0.945 on three held-out donors.

Each model was screened in four in-silico perturbation arms (rescue and
induction × delete and overexpress). Pathway drivers went to Enrichr
over-representation; a stricter high-confidence target set went to
connectivity-map drug-signature enrichment, intersected across the two models.

> **Seven scripts here are reconstructions.** The MT model's own pipeline was
> deleted from the cluster before this repository was assembled and no copy
> survives — nor does the MT checkpoint, so these scripts cannot currently be
> re-run even with the data in hand. They are the six files in
> `*/mt_multitask_RECONSTRUCTED/` plus `figures/mt_cls_umap_RECONSTRUCTED.py`.
> See the repository README and each file's header. All other scripts in this
> section are original.

## 00_environment/

Slurm setup jobs: conda environment, Geneformer package install,
Geneformer-V2-104M weight download, and a `transformers` dependency fix.

## 01_dataset_preparation/

Shared by both models.

| Script | What it does |
|---|---|
| `organize_10x_from_geo_tar.sh` | Unpacks the GEO tarball into per-GSM 10x triplet folders. |
| `parse_geo_series_metadata.py` | Parses `GSE209552_family.soft` into a per-sample TSV (donor, condition, region, age, sex). |
| `build_geneformer_h5ad_from_10x.py` | Merges the per-sample count matrices into one raw-counts `.h5ad`. |
| `qc_sweep.py` | Sweeps MAD/QC configurations against the published 15,110-nucleus target to find the one that reproduces the cohort's own filtering. |
| `rebuild_with_paper_qc.py`, `refine_qc_to_paper.py` | Rebuild under the selected QC. |
| `apply_ngene1000_filter.py`, `build_paper_matched_primary.py` | Apply the ≥1,000-gene cutoff. This alone reproduces the published QC set to within 62 nuclei (0.4%). Produces the 15,172-nucleus primary object. |
| `annotate_cell_types_de_novo.py` | Leiden clustering plus marker scoring over the canonical cortical marker panels; Scrublet doublet screening in the first pass. Nuclei admitted only by the final gene filter are labelled by repeating the procedure. |
| `create_geneformer_subsets.py`, `rebuild_subsets_and_tokenize.py` | Carve the cell-type-specific `.h5ad` files. |
| `tokenize_paperQC_nGene1000.py` | Rank-value encoding with `TranscriptomeTokenizer` (`model_version: V2`). |
| `validate_geneformer_h5ad.py` | Checks the `.h5ad` files against the tokenizer contract. |
| `01_count_neuronal_genes.py` | Gene-detection counts per neuronal class; the basis for the eligibility universe. |

## 02_model_finetuning/

**`ft_excitatory/`** — original code.

| Script | What it does |
|---|---|
| `12_finetune_tbi_classifier_excitatory.py` + `run_*.sh` | Fine-tunes the FT classifier. Donor-level test split, never seen in training; cell-level random splits would leak donor identity. |
| `12b_evaluate_finetuned_classifier.py` + `run_*.sh` | Eval-only pass on the held-out outer test set. |
| `13_extract_state_embs_finetuned.py` + `run_*.sh` | TBI and control state-embedding centroids from the FT model. |
| `finetuned_celltype_config.py` | Per-cell-type paths and labels shared by the ISP stages. |

**`mt_multitask_RECONSTRUCTED/`** — reconstructions, in run order:
`mt_01_depth_match_garza.py` (quantile-binned depth matching within cell type),
`mt_02_donor_split_garza.py` (one donor held out per condition; 80/20 train/valid
at seed 42), `mt_03_finetune_garza.py` (`MTLClassifier`, three heads, the
published hyperparameters, and the checkpoint-selection rule).

## 03_in_silico_perturbation/

Both models were run in the same four arms:

| Arm | Start cells | Goal | Perturbation | A positive hit means |
|---|---|---|---|---|
| Rescue / delete | TBI | control | delete | losing the gene moves TBI cells toward control |
| Rescue / overexpress | TBI | control | overexpress | activating the gene moves TBI cells toward control |
| Induction / delete | control | TBI | delete | losing the gene drives control cells toward TBI |
| Induction / overexpress | control | TBI | overexpress | activating the gene drives control cells toward TBI |

**`ft_excitatory/`** — original code. `perturb_cells_finetuned_*.py` with the
`run_perturb_array_*.sh` Slurm array wrappers run the perturbations; all
start-state nuclei are perturbed (2,047 TBI for rescue, 3,978 control for
induction), partitioned into fixed chunks of 56 nuclei to fit the wall-time
limit. `05_get_perturb_stats_finetuned.py` aggregates chunk output into per-gene
cosine-shift statistics. `perturb_safeguards.py` holds the statistical core —
per-gene aggregation, donor consistency, the 20%-detection eligibility rule,
coverage gates — and is imported by the MT aggregation too.
`07_biological_candidate_filter_finetuned.py` with `biological_gene_classes.py`
applies the biological-candidate filter.

**`mt_multitask_RECONSTRUCTED/`** — `mt_04_build_state_embs.py` (exact-mean CLS
centroids from the training split), `mt_05_in_silico_perturbation.py` (the four
arms), `mt_06_aggregate_isp.py` (per-gene statistics, stratified by cell type,
importing the surviving `perturb_safeguards.py`).

### A note on the MT screen's safeguards

The MT split holds out one donor per condition, so the 1,846 held-out TBI nuclei
all come from a single donor. Geneformer's donor-consistency score is therefore
undefined (NaN) for every gene in the MT screen and is neither filtered on nor
reported. The other safeguards are unaffected. Holding out a second control was
not possible: the cohort has three, and training on one would gut the disease
head.

MT state-embedding centroids come from the **training** split while perturbation
runs on **held-out** nuclei. This departs from Theodoris et al., who draw both
from the same pool. It avoids leakage between the centroid and the perturbed
cells, but it is a departure and is flagged as such in the Methods.

MT checkpoint selection also used the held-out donors rather than a validation
split, so the reported MT metrics are optimistically biased relative to the FT
metrics. The paper states this.

## 04_pathway_enrichment/

**Supplementary Figure 6.** `09_pathway_enrichment_finetuned.py` runs the
over-representation analysis (GSEApy `enrichr`, one-sided Fisher exact,
FDR-corrected within library) against GO Biological Process 2026, KEGG 2026,
Reactome 2024 and MSigDB Hallmark 2020. Background is the 7,374 genes passing
the 20%-detection filter in excitatory neurons, not the whole genome, so that
neuron-expressed pathways are not inflated.
`excit_pathway_heatmap_enrich.py` runs the unified Enrichr pass across arms;
`excit_pathway_heatmap_FT.py` and `excit_pathway_heatmap_MT.py` draw the two
heatmap panels.

Arms differ in foreground size by two orders of magnitude, so they differ in
term specificity. The FT induction/overexpression arm puts 2,621 foreground
genes against the 7,374-gene background and returns broad, weakly enriched terms;
it is read as descriptive. The rescue arms return few, specific terms.

## 05_drug_repurposing/

**Figure 2b.** `excit_drug_repurpose.py` assigns high-confidence targets to
direction-aware sets — INHIBIT (activation induces the TBI state, or deletion
rescues it) and ACTIVATE (the converse) — drops direction-ambiguous genes, and
tests INHIBIT targets against drug down-regulation signatures and ACTIVATE
targets against up-regulation signatures, following connectivity-map logic.
`excit_drug_crossmodel.py` intersects the two models' ranked tables and marks
the consensus compounds. `excit_drug_repurpose_figure.py` draws the panel.

Two scoring statistics are in play and are not the same: compounds are **ranked**
by the sum of −log10 FDR-corrected p across signatures, while Figure 2b
**plots** the best individual signature per model.

No LINCS L1000 signature passed FDR correction, so every candidate compound
derives from the GEO drug-perturbation libraries after exclusion of toxicants and
chemotherapy agents.

## 06_confounders_and_robustness/

| Script | What it does |
|---|---|
| `confound_check_sex_age_drivers.py` | Screens every pathway-driver list (four arms × two models) for Y-chromosome and X-inactivation genes and for genes whose excitatory-neuron expression tracks donor age. Age association is the Spearman correlation between donor age and pseudobulk log-CPM across the 10 TBI donors with ≥20 excitatory nuclei. |
| `cross_direction_network_analysis.py` + `run_cross_direction_network_v2.sh` | Relationship between the induction and rescue screens; the basis for the bidirectional core. |
| `compare_singlect_vs_jointpg.py` | Single-cell-type screen against the multi-task screen on the same cohort and the same arm. |

Disease state is confounded with donor covariates in this cohort — all three
controls are men aged 69–87, TBI donors are younger and include two women, and
each library comes from one donor, so batch is nested within donor. The analysis
contains these effects rather than removing them. Because the controls lie above
most of the TBI age range, age effects that coincide with disease state itself
cannot be separated.

## figures/

| Script | Output |
|---|---|
| `mt_cls_umap_RECONSTRUCTED.py` | **Figure 2a**. MT second-to-last-layer CLS embeddings of all 9,428 depth-matched nuclei, colored by disease state and outlined by cell type. **Reconstruction** — the original is among the deleted MT scripts. |
| `umap_all_celltypes_expression.py` | **Supplementary Figure 5a**. Gene-expression UMAP of all 15,172 QC-passing nuclei (2,000 HVGs, 30 PCs), colored by annotated cell type. Standard scanpy workflow on raw counts, no batch correction, so donor structure stays visible and judgeable. |
| `umap_excit_ft.py` | **Supplementary Figure 5b**. FT excitatory CLS embeddings by condition. |
| `make_suppfig5_model_metrics.py` | **Supplementary Figure 5**, metric panels. Confusion matrices and metric tables for both models. |
| `make_fig2_panelB.py` | **Figure 2b**. Cross-model consensus compounds. |
| `heatmap_excit_finetuned_cosine_shifts.py` | Per-cell cosine shifts for the top FT rescue drivers. |

Embeddings are the second-to-last-layer CLS token (768 dimensions), reduced with
UMAP (`n_neighbors` 30, `min_dist` 0.3, cosine metric, `random_state` 0).
Note that `emb_layer: -1` means the second-to-last layer in Geneformer's
convention, not the last; the classification heads read the final layer.

Source data for this section is in [`../data/02_geneformer_in_silico_screen`](../data/02_geneformer_in_silico_screen); see
[`../data/README.md`](../data/README.md) for formats and column schemas.
