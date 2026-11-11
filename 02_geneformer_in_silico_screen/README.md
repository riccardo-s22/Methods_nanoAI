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

> The seven `*_RECONSTRUCTED` scripts were rewritten from the Methods after the
> original MT code was lost; they document the procedure but were not used to
> generate the published results. All other scripts in this section are original.

## Rerunning this section

From a fresh clone, on an ordinary CPU machine, from any working directory:

    bash 02_geneformer_in_silico_screen/run_section.sh

This takes about one minute. It reads only `data/02_geneformer_in_silico_screen/`, never writes to
`data/`, and writes to `results/02_geneformer_in_silico_screen/` (one log per step in `logs/<step>.log`).
Stages 07 and 09 are part of the FT pipeline and write, as they always have,
into the pipeline working tree `work/Geneformer_TBI/results/` (`$GF_WORK`), where
the cross-direction analysis reads them. `results/` and `work/` are git-ignored.
Every step prints one `PASS`, `FAIL` or `SKIP` line, repeated in a summary at the
end; a failing step does not stop the run, and the exit status is non-zero if any
step failed. Three `SAME`/`DIFF` lines compare the regenerated source tables of
Supplementary Figure 6 and Figure 2b byte for byte with the copies in `data/`.
In the default run the two Supplementary Figure 6 matrices are read from `data/`
and written back by the heatmap scripts, so their `SAME` lines show only that
the scripts pass the shipped values through unchanged; they are not an
independent recomputation. The consensus table is recomputed from the shipped
per-model drug tables. Recomputing the matrices from the per-gene tables needs
`--with-project-tree` (see below).

| Option | Effect |
|---|---|
| `--offline` (or `GF_OFFLINE=1`) | Skip the five steps that query a web service (stage 09: Enrichr; cross-direction analysis: STRING). By default they run. |
| `--with-project-tree` (or `GF_WITH_PROJECT_TREE=1`) | Also run the class B steps below, when their inputs are found under `$GF_WORK`, `$GF_JOINTPG` and `$GF_MECHANO`. A step whose inputs are absent prints `SKIP` and names the first missing file. |
| `--only REGEX` | Run only the steps whose name matches, e.g. `--only heatmap`. Step names are the names in the PASS/FAIL/SKIP lines. |
| `--print-cluster-plan` | Print the upstream (cluster) half of the section as commands in run order and exit. Nothing is executed. |

### What can be rerun where

Every script of the section is listed in a table under its folder below, with a
class, the exact command (run from the repository root), what it reads, what it
writes and the time it took on the machine used for the port (8 cores, 9.9 GB
RAM, WSL2, no GPU).

| Class | Meaning | Run by `run_section.sh` |
|---|---|---|
| **A** | Runs from a fresh clone on a CPU, from the tables in `data/`. **A·net** also needs internet access. | yes |
| **B** | CPU-only, but needs files that are not shipped in `data/` (named in each table). | only with `--with-project-tree`, when the files exist |
| **C** | Needs a GPU, a model checkpoint or a Slurm cluster, or cannot run at all (the MT reconstructions). | never; see [Cluster walkthrough](#cluster-walkthrough) |

What the default run regenerates from `data/`: the two Supplementary Figure 6
heatmaps and their source matrices, the cross-model consensus table and both
renderings of Figure 2b, the Supplementary Figure 5 metric panels, the
biological-candidate tables (stage 07), their Enrichr tables (stage 09) and the
cross-direction analysis.

What it cannot regenerate: anything that starts from the snRNA-seq counts, a
model or the raw perturbation output. In particular the two steps that turn the
per-gene tables into the pathway-driver and drug-target gene sets
(`excit_pathway_heatmap_enrich.py`, `excit_drug_repurpose.py`) are class B: they
need the 20%-detection table that defines the 7,374-gene universe, which is not
shipped, and for the MT model the MT per-gene tables, which are not available.
The default run therefore starts from the shipped results of those two steps
(`pathway_enrichment/` and `drug_signature_enrichment/`). The 7,374-gene
universe, the pathway-driver lists and the INHIBIT/ACTIVATE target sets cannot
be re-derived from `data/` alone.

### Configuration

Nothing needs to be set for the default run. All locations come from
[`gf_paths.py`](gf_paths.py) (Python) and [`gf_env.sh`](gf_env.sh) (shell), which
resolve the same variables to the same defaults and document every one of them;
`python3 02_geneformer_in_silico_screen/gf_paths.py` prints the resolved configuration.

| Variable | Default | Meaning |
|---|---|---|
| `NANOAI_DATA` | `<repo>/data` | Shipped source data (read-only). |
| `NANOAI_RESULTS` | `<repo>/results` | Outputs of the CPU analyses and all logs. |
| `PYTHON` | `python3` | Interpreter of the class A steps. Needs numpy, pandas, scipy, matplotlib, requests and gseapy. |
| `GF_WORK` | `<repo>/work/Geneformer_TBI` | Working tree of the FT pipeline, in the layout of the original project (`h5ad/`, `metadata/`, `reports/`, `tokenized_paperQC_nGene1000/`, `models/Geneformer/`, `results/`). Pointing it at a copy of the original tree makes every script read what it originally read. |
| `GF_JOINTPG` | `<repo>/work/Geneformer_JointPG` | Tree of the MT model. The detection table and the MT per-gene tables are looked up below it. |
| `GF_MECHANO` | `<repo>/work/MechanoGEPred` | Only `data/mechano_genes_list.txt` is read from it (`compare_singlect_vs_jointpg.py`). |
| `GF_GEO` | `<repo>/work/GSE209552` | GEO download (stage 01). |
| `GF_MODEL_DIR` | `$GF_WORK/models/Geneformer` | Snapshot of `ctheodoris/Geneformer` (weights, token dictionaries). |
| `GF_FT_MODEL` | unset | A fine-tuned FT checkpoint to use directly, instead of the path in `$GF_WORK/results/finetune_tbi_excitatory/trained_model_path.txt`. |
| `GF_MT_MODEL` | `$GF_JOINTPG/models/garza_mt/ksplit1/checkpoint-600` | The retained MT checkpoint (not available). |
| `GF_CONDA_ENV`, `GF_PYTHON` | unset, `python` | Environment and interpreter with Geneformer and PyTorch: every Slurm script, and the class B/C steps that need `geneformer`, `pyarrow` or `umap-learn`. |
| `GF_PREP_PYTHON` | `$GF_PYTHON` | Interpreter with scanpy, anndata, harmonypy, scrublet and statsmodels: stage 01 and the class B steps that read an AnnData. |
| `STRING_CALLER_IDENTITY` | `geneformer_tbi_in_silico_screen` | Client name sent to the STRING API by the cross-direction analysis. |

`gf_env.sh` lists the remaining variables (`GF_SOFT`, `GF_SRA_RUN_TABLE`,
`GF_PRETRAINED`, `GF_GENEFORMER_PKG`, `GF_FT_MODELS`, `GF_LOGS`).

Where a CPU analysis reads the output of an earlier one, it looks in three
places, in this order: `results/02_geneformer_in_silico_screen/` (a rerun in this repository), then
`$GF_WORK/results/` (the original project tree), then the copy shipped in
`data/`. A rerun therefore feeds the next step, and a fresh clone falls back on
the shipped tables.

### Checks made when the section was ported

Reruns on the port machine (Python 3.13.13, NumPy 2.5.1, pandas 3.0.3,
matplotlib 3.11.0, GSEApy 1.3.1; the Methods give GSEApy 1.2.1), compared with
the shipped tables and with the output files of the original run:

- **Supplementary Figure 6 matrices** (`excit_pathway_heatmap_FT.py`, `_MT.py`):
  byte-identical to `data/.../pathway_enrichment/`. The FT matrix holds 7,837
  terms, 367 with FDR-corrected p < 0.05 in at least one arm (2, 7, 6 and 354 per
  arm); the MT matrix holds 37 terms (2, 17 and 18 in its three tested arms).
  The matrices hold FDR-corrected p values; the scripts take −log10 when drawing.
- **Cross-model consensus** (`excit_drug_crossmodel.py`): byte-identical to
  `data/.../crossmodel_consensus.tsv` (eight rows; six compounds after `vx` and
  `sodium` are dropped by the figure scripts).
- **Stage 07**: the eight candidate tables are byte-identical to the original
  run's.
- **Stage 09** (Enrichr, queried at the time of the port): all 38,369 term rows of the 32 library
  tables have the same p values, FDR-corrected p values, odds ratios and gene
  lists as the original run. Rows with tied FDR-corrected p are written in a
  different order, so 26 of the 32 files and two of the four `summary_top5.tsv`
  are not byte-identical.
- **Cross-direction analysis**: the overlap tables (`01`, `02`) are
  byte-identical and the pathway-convergence table (`03`) has the same content;
  the STRING tables (`04`–`08`) differ, because STRING now returns 485 edges
  where the original run got 392 (326 in common).
- **Figures**: the two heatmaps and the three-panel drug figure were compared
  with the original PNG files by eye, not byte for byte: the regenerated images
  are 1 to 4 pixels larger (matplotlib 3.11.0 here, 3.10.5 in the Methods). The
  heatmaps drawn from the shipped matrices are byte-identical to those drawn
  from the original per-arm Enrichr tables on the same machine. The two
  `make_*` scripts ask for Arial and fall back to DejaVu Sans where it is not
  installed, as on the port machine.

Class B steps, rerun against a copy of the original project tree:

- `build_detection_table.py` gives 7,374 genes detected in ≥20% of either group.
  With that table, `excit_pathway_heatmap_enrich.py` returns FT pathway-driver
  lists byte-identical to the original run's (2,192, 1,303, 2,410 and 2,621
  genes) and Enrichr statistics identical for every term. The regenerated FT
  matrix differs from the shipped one in three Reactome term names only: the
  shipped file has `IKK?` where Enrichr returns `IKKε`.
- `excit_drug_repurpose.py FT union|core`: target sets of 1,809 INHIBIT and
  1,471 ACTIVATE genes, bidirectional core 594 and 596, 14 direction-conflict
  genes dropped; Enrichr statistics identical for all terms. The core table and
  `all_sig_drug_terms.tsv` are byte-identical to the shipped ones; the union
  table differs in the sixteenth significant digit of `sum_signal` for 9 of its
  45 rows.
- `confound_check_sex_age_drivers.py`: `age_corr.tsv` is byte-identical to the
  original; `confound_summary.tsv` is identical except that the original file
  also lists `KDM5C` in two rows, a gene the script in this repository does not
  search for.
- `umap_all_celltypes_expression.py`: coordinates identical to the original
  file for all 15,172 nuclei. `umap_excit_ft.py`, run from the original
  embedding cache: 6,025 nuclei (2,047 TBI, 3,978 control); the UMAP coordinates
  are not numerically identical to the original ones (umap-learn 0.5.12 here,
  0.5.9.post2 in the Methods).
- `heatmap_excit_finetuned_cosine_shifts.py`, run on the original chunk pickles
  and the rebuilt detection table: the 40 x 2,047 shift matrix is identical to
  the original one.

The MT half of these steps, `compare_singlect_vs_jointpg.py` and the seven
`*_RECONSTRUCTED` scripts could not be rerun: the MT tables and checkpoint are
not available.

### For coding agents

- Run `bash 02_geneformer_in_silico_screen/run_section.sh` and read the summary; a per-step log is in
  `results/02_geneformer_in_silico_screen/logs/`. Use `--offline` when there is no network.
- Never write under `data/`. Do not change a statistic, threshold, seed or
  library name to make an output match; report the difference instead.
- Do not try to run class C steps on a machine without a GPU. Print them with
  `--print-cluster-plan`.
- A `DIFF` line is expected only after inputs were regenerated (for example
  with `--with-project-tree`); on a fresh clone all three lines read `SAME`.

## 00_environment/

Slurm setup jobs: conda environment, Geneformer package install,
Geneformer-V2-104M weight download, and a `transformers` dependency fix.

**Run.** Slurm jobs; submit from the repository root after `mkdir -p results/02_geneformer_in_silico_screen/logs`. They need network access and an existing Python 3.10 environment with CUDA PyTorch, given as `GF_CONDA_ENV` (or `GF_PYTHON`). No script creates that environment; its recipe is [`../envs/geneformer.yml`](../envs/geneformer.yml), whose pins are transcribed from the paper and were not tested (read its header). Without Slurm: `SLURM_JOB_ID=0 SLURMD_NODENAME=local bash <script>`. Not rerun for the port.

| Script | Class | Command | Reads | Writes | Time / resources |
|---|---|---|---|---|---|
| `00_install_and_download.sh` | C | `sbatch 02_geneformer_in_silico_screen/00_environment/00_install_and_download.sh` | PyPI | packages in the environment | 1st of 4; CPU partition, 4 CPUs, 16–32 GB, 1–3 h limits |
| `00c_install_geneformer_pkg.sh` | C | `sbatch 02_geneformer_in_silico_screen/00_environment/00c_install_geneformer_pkg.sh` | Hugging Face `ctheodoris/Geneformer` (clone without LFS) | `geneformer` package in the environment | 2nd |
| `00e_fix_transformers_deps.sh` | C | `sbatch 02_geneformer_in_silico_screen/00_environment/00e_fix_transformers_deps.sh` | PyPI | `transformers` 4.46.* pinned | 3rd |
| `00d_download_model_weights.sh` | C | `sbatch 02_geneformer_in_silico_screen/00_environment/00d_download_model_weights.sh` | Hugging Face | `$GF_MODEL_DIR` (Geneformer-V2-104M, token dictionaries; 2.8 GB) | 4th |
| `00_install_geneformer_env.sh`, `00b_download_geneformer_model.sh` | C | not needed |  |  | failed in the original run and were superseded by the four above; kept as the record |

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

**Run.** CPU, in the order of the table; inputs come from GEO, not from `data/`. The scanpy scripts use `$GF_PREP_PYTHON`, the tokenizer `$GF_PYTHON`. Three things are not covered by an original script: the GEO download (now `fetch_geo_GSE209552.sh`, added for this repository), the creation of the two conda environments, and the copy of each cell-type subset into `h5ad/geneformer_inputs_paperQC_nGene1000/<type>/`, which the tokenizer reads (the command is printed by `run_section.sh --print-cluster-plan`). `parse_geo_series_metadata.py` does not regenerate the `region` value of the 12 TBI samples (it leaves it empty), so a rebuilt object lacks `region` for the TBI nuclei; counts and every other column are unaffected, and the FT model does not use `region`.

| Script | Class | Command | Reads | Writes | Time / resources |
|---|---|---|---|---|---|
| `fetch_geo_GSE209552.sh` | B (network) | `bash 02_geneformer_in_silico_screen/01_dataset_preparation/fetch_geo_GSE209552.sh` | GEO GSE209552 | `$GF_GEO/GSE209552_RAW/`, `$GF_GEO/GSE209552_family.soft`, `$GF_WORK/reports/snrna_sample_manifest.txt` | 42 s; 0.9 GB download. Not part of the original run |
| `parse_geo_series_metadata.py` | B | `$GF_PREP_PYTHON 02_geneformer_in_silico_screen/01_dataset_preparation/parse_geo_series_metadata.py` | `$GF_SOFT`; optional `$GF_SRA_RUN_TABLE` | `$GF_WORK/metadata/gse209552_sample_metadata_curated.tsv` | 3 s |
| `build_geneformer_h5ad_from_10x.py` | B | `$GF_PREP_PYTHON 02_geneformer_in_silico_screen/01_dataset_preparation/build_geneformer_h5ad_from_10x.py` | `$GF_GEO/GSE209552_RAW/`, the manifest, the metadata table | `$GF_WORK/h5ad/Garza_GSE209552_all_raw_counts_geneformer_ready.h5ad` (24,621 nuclei) | 1 min 49 s; 1.8 GB RAM |
| `rebuild_with_paper_qc.py` | B | `$GF_PREP_PYTHON 02_geneformer_in_silico_screen/01_dataset_preparation/rebuild_with_paper_qc.py` | the merged `.h5ad` | `$GF_WORK/h5ad/..._paperQC_annotated.h5ad` | not rerun |
| `qc_sweep.py` | B | `$GF_PREP_PYTHON 02_geneformer_in_silico_screen/01_dataset_preparation/qc_sweep.py` | the two `.h5ad` above | `$GF_WORK/reports/qc_sweep.tsv` | 56 s |
| `build_paper_matched_primary.py` | B | `$GF_PREP_PYTHON 02_geneformer_in_silico_screen/01_dataset_preparation/build_paper_matched_primary.py` | the two `.h5ad` above | `$GF_WORK/h5ad/..._paperQC_nGene1000_annotated.h5ad` (15,172 nuclei), `h5ad/celltype_subsets_paperQC_nGene1000/`, `reports/paperQC_nGene1000_report.txt` | needs more than 8 GB RAM (stopped at 8.1 GB on the port machine) |
| `tokenize_paperQC_nGene1000.py` | B | `$GF_PYTHON 02_geneformer_in_silico_screen/01_dataset_preparation/tokenize_paperQC_nGene1000.py` | `$GF_WORK/h5ad/geneformer_inputs_paperQC_nGene1000/<type>/` | `$GF_WORK/tokenized_paperQC_nGene1000/<type>/*.dataset` | not rerun |
| `organize_10x_from_geo_tar.sh`, `annotate_cell_types_de_novo.py`, `create_geneformer_subsets.py`, `rebuild_subsets_and_tokenize.py`, `validate_geneformer_h5ad.py`, `01_count_neuronal_genes.py`, `apply_ngene1000_filter.py`, `refine_qc_to_paper.py` | B | not part of the final sequence |  |  | earlier attempts kept as the record. The last two write the same output names as `build_paper_matched_primary.py`: do not run them after it |

## 02_model_finetuning/

**`ft_excitatory/`** — original code.

| Script | What it does |
|---|---|
| `12_finetune_tbi_classifier_excitatory.py` + `run_*.sh` | Fine-tunes the FT classifier. Donor-level test split, never seen in training; cell-level random splits would leak donor identity. |
| `12b_evaluate_finetuned_classifier.py` + `run_*.sh` | Eval-only pass on the held-out outer test set. |
| `13_extract_state_embs_finetuned.py` + `run_*.sh` | TBI and control state-embedding centroids from the FT model. |
| `finetuned_celltype_config.py` | Per-cell-type paths and labels shared by the ISP stages. |

**`mt_multitask_RECONSTRUCTED/`** — in run order:
`mt_01_depth_match_garza.py` (quantile-binned depth matching within cell type),
`mt_02_donor_split_garza.py` (one donor held out per condition; 80/20 train/valid
at seed 42), `mt_03_finetune_garza.py` (`MTLClassifier`, three heads, the
published hyperparameters, and the checkpoint-selection rule).

**Run.** `ft_excitatory/` needs one GPU with more than 11 GB of memory, the tokenized excitatory dataset and the pretrained weights; submit the `run_*.sh` scripts from the repository root. Not rerun for the port. The `mt_multitask_RECONSTRUCTED/` scripts were not used for the published results and cannot run end to end; their default locations (`gf_paths.MT_RECON`) are a layout chosen for this repository, not that of the original MT project.

| Script | Class | Command | Reads | Writes | Time / resources |
|---|---|---|---|---|---|
| `12_finetune_tbi_classifier_excitatory.py` | C | `sbatch 02_geneformer_in_silico_screen/02_model_finetuning/ft_excitatory/run_finetune_tbi_classifier_excitatory.sh` | tokenized excitatory dataset, `$GF_PRETRAINED` | `$GF_WORK/results/finetune_tbi_excitatory/` (incl. `heldout_eval/`, `trained_model_path.txt`), `$GF_FT_MODELS/<date>_geneformer_cellClassifier_excit_tbi_clf/` | 1 GPU, 8 CPUs, 64 GB, 6 h limit; original run 306 steps in 995 s |
| `12b_evaluate_finetuned_classifier.py` | C | `sbatch 02_geneformer_in_silico_screen/02_model_finetuning/ft_excitatory/run_evaluate_finetuned_classifier.sh` | an existing checkpoint (`GF_FT_MODEL`), the labelled test set | `$GF_WORK/results/finetune_tbi_excitatory/heldout_eval/` | 1 GPU, 4 CPUs, 32 GB, 30 min limit |
| `13_extract_state_embs_finetuned.py` | C | `sbatch 02_geneformer_in_silico_screen/02_model_finetuning/ft_excitatory/run_extract_state_embs_finetuned.sh` | checkpoint, tokenized excitatory dataset | `$GF_WORK/results/perturb_neurons_finetuned_control_to_tbi/state_embs/excitatory_finetuned_state_embs.pkl` | 1 GPU, 4 CPUs, 32 GB; original run 5.4 min |
| `finetuned_celltype_config.py` | library | imported by stages 05, 07 and 09 |  |  |  |
| `mt_01_depth_match_garza.py` | B | `$GF_PREP_PYTHON 02_geneformer_in_silico_screen/02_model_finetuning/mt_multitask_RECONSTRUCTED/mt_01_depth_match_garza.py` | the 15,172-nucleus `.h5ad` | `$GF_JOINTPG/h5ad/garza_mt_depth_matched.h5ad` | reconstruction; not run |
| `mt_02_donor_split_garza.py`, `mt_03_finetune_garza.py` | C | see the script headers | depth-matched corpus; tokenized MT corpora that no script writes | `$GF_JOINTPG/...` | reconstruction; cannot run end to end |

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

**Run.** The perturbation itself is 218 GPU array tasks (about 865 GPU-hours in the original logs) and writes 1.4 GB of chunk pickles; stage 05 aggregates them on a CPU. Neither can be rerun from `data/`, which ships the output of stage 05. Stage 07 can: `run_section.sh` first copies the shipped per-gene tables into the working-tree layout (`python3 02_geneformer_in_silico_screen/gf_paths.py --seed-ft-stats`, which never overwrites a table already there) and then runs it for the four arms.

| Script | Class | Command | Reads | Writes | Time / resources |
|---|---|---|---|---|---|
| `perturb_cells_finetuned_tbi_to_control.py` | C | `sbatch 02_geneformer_in_silico_screen/03_in_silico_perturbation/ft_excitatory/run_perturb_array_finetuned_tbi_to_control_delete.sh` and `..._overexpress.sh` | checkpoint, state embeddings, tokenized excitatory dataset | `$GF_WORK/results/perturb_neurons_finetuned_tbi_to_control/chunks_<mode>/` | 37 array tasks per mode; 1 GPU, 8 CPUs, 32 GB, 12 h limit each |
| `perturb_cells_finetuned_control_to_tbi.py` | C | `sbatch 02_geneformer_in_silico_screen/03_in_silico_perturbation/ft_excitatory/run_perturb_array_finetuned_control_to_tbi_delete.sh` and `..._overexpress.sh` | as above | `$GF_WORK/results/perturb_neurons_finetuned_control_to_tbi/chunks_<mode>/` | 72 array tasks per mode; same resources |
| `05_get_perturb_stats_finetuned.py` | B | `$GF_PYTHON 02_geneformer_in_silico_screen/03_in_silico_perturbation/ft_excitatory/05_get_perturb_stats_finetuned.py --mode <delete|overexpress> --celltype <excit_rescue|excit>`; Slurm: `run_05_get_perturb_stats_finetuned_delete.sh`, `..._overexpress.sh` (induction), `run_05_07_09_finetuned_rescue_overexpress.sh` | chunk pickles, tokenized excitatory dataset, token dictionaries (none shipped) | `$GF_WORK/results/perturb_neurons_finetuned_<direction>/stats_<mode>/*_unfiltered_all_gene_ranked.tsv`, `*_filtered_high_confidence_gene_ranked.tsv` (the tables in `data/`) | 41 min and 1.5 GB RAM for one arm on the port machine; 13–14 min on the cluster (8 CPUs, 64 GB) |
| `07_biological_candidate_filter_finetuned.py` | A | `python3 02_geneformer_in_silico_screen/gf_paths.py --seed-ft-stats`, then `python3 02_geneformer_in_silico_screen/03_in_silico_perturbation/ft_excitatory/07_biological_candidate_filter_finetuned.py --mode <delete|overexpress> --celltype <excit_rescue|excit>` | `data/{S}/in_silico_perturbation/ft_excitatory/<arm>/*_unfiltered_all_gene_ranked.tsv` (through the seeded copy) | `$GF_WORK/results/perturb_neurons_finetuned_<direction>/stats_<mode>/*_filtered_biological_candidate_ranked.tsv`, `*_top100_biological_candidate_summary.tsv` | 1–2 s per arm |
| `perturb_safeguards.py`, `biological_gene_classes.py` | library | imported |  |  |  |
| `mt_04_build_state_embs.py`, `mt_05_in_silico_perturbation.py`, `mt_06_aggregate_isp.py` | C | see the script headers | MT checkpoint and tokenized MT corpora (not available) | `$GF_JOINTPG/results/isp_RECONSTRUCTED/` | reconstruction; cannot run |

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

**Run.** Two different enrichment passes live in this folder. `excit_pathway_heatmap_enrich.py` is the one behind Supplementary Figure 6: Enrichr libraries `GO_Biological_Process_2026`, `KEGG_2026`, `Reactome_Pathways_2024` and `MSigDB_Hallmark_2020`, pathway-driver genes as foreground, the 7,374-gene universe as background. `09_pathway_enrichment_finetuned.py` is the earlier pipeline stage that follows stage 07: as written it queries `GO_Biological_Process_2023`, `Reactome_2022`, `KEGG_2021_Human` and `MSigDB_Hallmark_2020` with the biological candidates as foreground and all tested genes of the arm as background; its output feeds only the cross-direction analysis.

The heatmap scripts read one `<arm>_enrichr_allterms.csv` per arm. When that file is found neither in a rerun nor in `$GF_WORK`, they take the arm's FDR-corrected p values from the shipped matrix, which the same script wrote from those files in the published run; the drawn heatmap is byte-identical either way.

| Script | Class | Command | Reads | Writes | Time / resources |
|---|---|---|---|---|---|
| `09_pathway_enrichment_finetuned.py` | A·net | `python3 02_geneformer_in_silico_screen/04_pathway_enrichment/09_pathway_enrichment_finetuned.py --mode <delete|overexpress> --celltype <excit_rescue|excit>` | stage 07 candidate table and the all-gene table of the arm, below `$GF_WORK/results/` | `$GF_WORK/results/perturb_neurons_finetuned_<direction>/enrichment_<mode>/{high_tier,all_bio}/<library>.tsv`, `summary_top5.tsv` | 8–9 s per arm |
| `build_detection_table.py` | B | `$GF_PREP_PYTHON 02_geneformer_in_silico_screen/04_pathway_enrichment/build_detection_table.py` | `$GF_WORK/h5ad/celltype_subsets_paperQC_nGene1000/excitatory_neurons_raw_geneformer_ready.h5ad` | `$GF_JOINTPG/results/isp/garza_disease_binary_TBI_to_control/mechano_analysis/layer3/garza_qc4_de_by_celltype.tsv` | 6 s. Not part of the original run: the original table came from the lost MT project |
| `excit_pathway_heatmap_enrich.py` | B·net | `python3 02_geneformer_in_silico_screen/04_pathway_enrichment/excit_pathway_heatmap_enrich.py [ARM ...]` (arms: `FT_delete FT_OE FT_del_c2t FT_oe_c2t MT_delete MT_OE MT_del_c2t MT_oe_c2t`; default all) | the detection table (not shipped); FT per-gene tables (shipped); MT per-gene tables (not available) | `results/02_geneformer_in_silico_screen/pathway_heatmap/enrich/<arm>/<arm>_hits_fdr05.tsv`, `<arm>_enrichr_allterms.csv` | 7 s for the four FT arms |
| `excit_pathway_heatmap_FT.py` | A | `python3 02_geneformer_in_silico_screen/04_pathway_enrichment/excit_pathway_heatmap_FT.py` | `data/02_geneformer_in_silico_screen/pathway_enrichment/excit_pathway_enrichment_matrix_FT.tsv` (or per-arm Enrichr tables) | `results/02_geneformer_in_silico_screen/pathway_heatmap/excit_pathway_enrichment_heatmap_FT.png`, `excit_pathway_enrichment_matrix_FT.tsv` | 2 s |
| `excit_pathway_heatmap_MT.py` | A | `python3 02_geneformer_in_silico_screen/04_pathway_enrichment/excit_pathway_heatmap_MT.py` | `data/02_geneformer_in_silico_screen/pathway_enrichment/excit_pathway_enrichment_matrix_MT.tsv` (or per-arm Enrichr tables) | `results/02_geneformer_in_silico_screen/pathway_heatmap/excit_pathway_enrichment_heatmap_MT.png`, `excit_pathway_enrichment_matrix_MT.tsv` | 3 s |

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

**Run.** `excit_drug_repurpose.py` takes the model and the mode as arguments and is run four times, union before core for each model: both modes write `<model>/all_sig_drug_terms.tsv`, and the shipped copy is the one left by the core run.

| Script | Class | Command | Reads | Writes | Time / resources |
|---|---|---|---|---|---|
| `excit_drug_repurpose.py` | B·net | `python3 02_geneformer_in_silico_screen/05_drug_repurposing/excit_drug_repurpose.py <FT|MT> <union|core>` | the detection table (not shipped); FT per-gene tables (shipped); MT per-gene tables (not available) | `results/02_geneformer_in_silico_screen/drug_repurpose/<model>/drug_candidates_ranked_<mode>.tsv`, `all_sig_drug_terms.tsv`, `<set>_<mode>/<set>_<mode>_allterms.csv` | 6–7 s per FT run |
| `excit_drug_crossmodel.py` | A | `python3 02_geneformer_in_silico_screen/05_drug_repurposing/excit_drug_crossmodel.py` | `data/02_geneformer_in_silico_screen/drug_signature_enrichment/{FT,MT}/drug_candidates_ranked_union.tsv`, `FT/drug_candidates_ranked_core.tsv` | `results/02_geneformer_in_silico_screen/drug_repurpose/crossmodel_consensus.tsv` | 1 s |
| `excit_drug_repurpose_figure.py` | A | `python3 02_geneformer_in_silico_screen/05_drug_repurposing/excit_drug_repurpose_figure.py` | the consensus table (the rerun above, else the shipped one), `FT/drug_candidates_ranked_core.tsv`, `FT/all_sig_drug_terms.tsv` | `results/02_geneformer_in_silico_screen/drug_repurpose/excit_drug_repurpose_figure.{png,pdf}` | 3 s |

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

**Run.**

| Script | Class | Command | Reads | Writes | Time / resources |
|---|---|---|---|---|---|
| `confound_check_sex_age_drivers.py` | B | `$GF_PREP_PYTHON 02_geneformer_in_silico_screen/06_confounders_and_robustness/confound_check_sex_age_drivers.py` | `$GF_WORK/h5ad/celltype_subsets_paperQC_nGene1000/excitatory_neurons_raw_geneformer_ready.h5ad` and the `<arm>_hits_fdr05.tsv` lists of `excit_pathway_heatmap_enrich.py` (neither shipped) | `results/02_geneformer_in_silico_screen/confounders/confound_summary.tsv`, `age_corr.tsv` | 8–17 s; 1.3 GB RAM. Needs anndata and statsmodels |
| `cross_direction_network_analysis.py` | A·net | `python3 02_geneformer_in_silico_screen/06_confounders_and_robustness/cross_direction_network_analysis.py` (Slurm: `run_cross_direction_network_v2.sh`) | stage 07 candidate tables, all-gene tables and stage 09 `enrichment_<mode>/` folders below `$GF_WORK/results/` | `results/02_geneformer_in_silico_screen/cross_direction_network_v2/01_...08_*.tsv` | 5 s. Queries STRING |
| `compare_singlect_vs_jointpg.py` | B·net | `python3 02_geneformer_in_silico_screen/06_confounders_and_robustness/compare_singlect_vs_jointpg.py` | zero-shot single-cell-type tables `$GF_WORK/results/perturb_<oligo|astro|microglia>_paperQC_nGene1000/stats/`, MT per-cell-type tables, `$GF_MECHANO/data/mechano_genes_list.txt` (none shipped; no script in the repository writes the first two) | `results/02_geneformer_in_silico_screen/confounders/singlect_vs_jointpg_comparison.tsv` | not rerun |

## figures/

| Script | Output |
|---|---|
| `mt_cls_umap_RECONSTRUCTED.py` | **Figure 2a**. MT second-to-last-layer CLS embeddings of all 9,428 depth-matched nuclei, colored by disease state and outlined by cell type. Rewritten; see note above. |
| `umap_all_celltypes_expression.py` | **Supplementary Figure 5a**. Gene-expression UMAP of all 15,172 QC-passing nuclei (2,000 HVGs, 30 PCs), colored by annotated cell type. Standard scanpy workflow on raw counts, no batch correction, so donor structure stays visible and judgeable. |
| `umap_excit_ft.py` | **Supplementary Figure 5b**. FT excitatory CLS embeddings by condition. |
| `make_suppfig5_model_metrics.py` | **Supplementary Figure 5**, metric panels. Confusion matrices and metric tables for both models. |
| `make_fig2_panelB.py` | **Figure 2b**. Cross-model consensus compounds. |
| `heatmap_excit_finetuned_cosine_shifts.py` | Per-cell cosine shifts for the top FT rescue drivers. |

Embeddings are the second-to-last-layer CLS token (768 dimensions), reduced with
UMAP (`n_neighbors` 30, `min_dist` 0.3, cosine metric, `random_state` 0).
Note that `emb_layer: -1` means the second-to-last layer in Geneformer's
convention, not the last; the classification heads read the final layer.

**Run.** Figures are written to `results/02_geneformer_in_silico_screen/figures/` and `results/02_geneformer_in_silico_screen/umap/`, never next to the scripts.

| Script | Class | Command | Reads | Writes | Time / resources |
|---|---|---|---|---|---|
| `mt_cls_umap_RECONSTRUCTED.py` | C | `$GF_PYTHON 02_geneformer_in_silico_screen/figures/mt_cls_umap_RECONSTRUCTED.py` | MT checkpoint and tokenized MT corpus (not available) | `$GF_JOINTPG/results/isp_RECONSTRUCTED/umap/` | reconstruction; cannot run |
| `umap_all_celltypes_expression.py` | B | `$GF_PREP_PYTHON 02_geneformer_in_silico_screen/figures/umap_all_celltypes_expression.py` | `$GF_WORK/h5ad/Garza_GSE209552_all_raw_counts_paperQC_nGene1000_annotated.h5ad` (not shipped) | `results/02_geneformer_in_silico_screen/umap/all_celltypes_expr_umap.{png,pdf}`, `all_celltypes_expr_umap_coords.tsv.gz` | 2 min; 2.3 GB RAM |
| `umap_excit_ft.py` | C (B with the embedding cache) | `$GF_PYTHON 02_geneformer_in_silico_screen/figures/umap_excit_ft.py` | FT checkpoint and tokenized excitatory dataset on a GPU; or an existing `excit_ft_cellembs.parquet` in the output folder or in `$GF_WORK/results/umap/` (none shipped) | `results/02_geneformer_in_silico_screen/umap/excit_ft_umap.png`, `excit_ft_umap_coords.parquet`, `excit_ft_cellembs.parquet` | about 1 min from the cache; 0.7 GB RAM. Needs pyarrow and umap-learn |
| `make_suppfig5_model_metrics.py` | A | `python3 02_geneformer_in_silico_screen/figures/make_suppfig5_model_metrics.py` | none: the metrics are written in the script, with their provenance | `results/02_geneformer_in_silico_screen/figures/suppfig5_metrics__*.{svg,png}` | 2 s |
| `make_fig2_panelB.py` | A | `python3 02_geneformer_in_silico_screen/figures/make_fig2_panelB.py` | the consensus table (the rerun of `excit_drug_crossmodel.py`, else `data/02_geneformer_in_silico_screen/drug_signature_enrichment/crossmodel_consensus.tsv`) | `results/02_geneformer_in_silico_screen/figures/Figure2b_crossmodel_consensus.{svg,png}` | 1 s |
| `heatmap_excit_finetuned_cosine_shifts.py` | B | `$GF_PYTHON 02_geneformer_in_silico_screen/figures/heatmap_excit_finetuned_cosine_shifts.py` | rescue/delete chunk pickles `$GF_WORK/results/perturb_neurons_finetuned_tbi_to_control/chunks_delete/`, the detection table, token dictionaries (none shipped) | `results/02_geneformer_in_silico_screen/figures/excit_finetuned_cosine_shift_heatmap.png`, `excit_finetuned_cosine_shift_matrix.npz` | 3 min 10 s reading the 4,801 pickles from a slow disk |

## Cluster walkthrough

The upstream half (stages 00 to 03) was run on a Slurm cluster with GPUs and is
not run by `run_section.sh`. The full sequence of commands, in run order, is
printed by

    bash 02_geneformer_in_silico_screen/run_section.sh --print-cluster-plan

In outline:

1. **Configure.** Export `GF_CONDA_ENV` (or `GF_PYTHON`) and `GF_PREP_PYTHON`,
   optionally `GF_WORK`, then `source 02_geneformer_in_silico_screen/gf_env.sh`. Slurm
   scripts are submitted from the repository root after
   `mkdir -p results/02_geneformer_in_silico_screen/logs`; they find the repository through
   `SLURM_SUBMIT_DIR` (or `NANOAI_REPO`, if exported), take the interpreter from
   `GF_PYTHON`, and write their logs to that folder. `sbatch` forwards the
   exported variables. No script carries a `--mail-user` line: pass your own
   with `sbatch --mail-user=<address>`. Partitions, constraints, node
   exclusions and limits in the `#SBATCH` headers are those of the original
   cluster.
2. **Stage 00, environment** (CPU, network): `00_install_and_download.sh`,
   `00c_install_geneformer_pkg.sh`, `00e_fix_transformers_deps.sh`,
   `00d_download_model_weights.sh`.
3. **Stage 01, dataset** (CPU, about 16 GB RAM): `fetch_geo_GSE209552.sh`,
   `parse_geo_series_metadata.py`, `build_geneformer_h5ad_from_10x.py`,
   `rebuild_with_paper_qc.py`, `qc_sweep.py`, `build_paper_matched_primary.py`,
   the manual copy into `geneformer_inputs_paperQC_nGene1000/`,
   `tokenize_paperQC_nGene1000.py`.
4. **Stage 02, FT model** (1 GPU): `run_finetune_tbi_classifier_excitatory.sh`,
   then `run_extract_state_embs_finetuned.sh`.
   `run_evaluate_finetuned_classifier.sh` is needed only to evaluate an existing
   checkpoint.
5. **Stage 03, perturbation** (GPU arrays): the four
   `run_perturb_array_finetuned_*.sh` scripts (rescue: 37 tasks per mode with
   `--constraint=a100-80g`; induction: 72 tasks per mode with
   `--constraint=bf16`), resubmitting tasks that hit the wall-time limit with
   `sbatch --array=<ids>`; then stage 05 for the four arms.
6. **CPU half.** `bash 02_geneformer_in_silico_screen/run_section.sh --with-project-tree`
   rebuilds the detection table, selects the pathway drivers and drug targets,
   queries Enrichr, and redraws the figures from the tree just built.

The MT model cannot be retrained from this repository: its checkpoint and
per-gene tables are not available, no script writes its tokenized corpora, and
the `*_RECONSTRUCTED` scripts document the procedure only.

Source data for this section is in [`../data/02_geneformer_in_silico_screen`](../data/02_geneformer_in_silico_screen); see
[`../data/README.md`](../data/README.md) for formats and column schemas.
