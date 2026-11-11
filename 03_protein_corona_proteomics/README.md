# Results 3 — the SWCNT protein corona

**Figure 3.**

The protein corona that forms on the (GT)15-SWCNTs is what determines each
condition's spectral features, so it was characterized directly. SWCNT–corona
complexes recovered after 24 h in blast-conditioned (n = 3) or sham-conditioned
(n = 3) medium were washed, the bound proteins eluted, and the eluate analysed
by directDIA in Spectronaut 17 against the human proteome (UniProt UP000005640).
6,198 protein groups were identified and 5,850 quantified in all six samples.

## ms_analysis.py

One script covers the whole analysis. Input is the Spectronaut PG report
(`../data/03_protein_corona_proteomics/Riccado_Report.csv`); outputs are the figure panels plus
`TBI_blast_proteomics_results.xlsx` and `summary.json`, written to
`../results/03_protein_corona_proteomics/`. Commands are in
[Running this section](#running-this-section).

Stages, in order:

1. **Canonical protein group selection** (`pick_canonical`) and UniProt
   accession → gene symbol mapping via the g:Profiler g:Convert REST API.
2. **Normalization.** log2 of the raw intensities, then median-centring across
   samples to a common median.
3. **Filtering.** A protein is kept if it is quantified in at least 2 of 3
   controls **and** at least 2 of 3 blast samples.
4. **Differential abundance** (`limma_fit`). Per-protein weighted least squares
   with empirical-Bayes moderation of the variances, as `limma`'s `lmFit` plus
   `eBayes`. **Missing values are handled natively, by fitting each protein on
   the samples where it was measured — they are not imputed.** A differentially
   abundant protein is one with BH-adjusted p < 0.05; no fold-change cutoff is
   applied. A batch-matched model and a no-batch-term model are fitted alongside
   as sensitivity checks.
5. **PCA** on the complete-case matrix of the 5,850 proteins quantified in every
   sample. Blast and control samples separate along PC1 (37.9% of variance).
6. **Functional enrichment** with g:Profiler (g:GOSt) through its REST API —
   separate from the Enrichr-based analysis of section 2 — run on the top 10%
   most abundant proteins against the tested proteome as background, and
   separately on the up- and down-regulated sets.
7. **Pre-ranked GSEA** on the moderated *t* statistic.
8. **Panels**: detected-proteome summary, PCA, volcano, top-protein heatmap and
   enrichment dot plots.

### Sample assignment

Control: C4R, C5R, C9. Blast: A1, A2, C6R. **C4 and C5 are excluded** — they
belong to a separate batch; C4R and C5R are their re-runs.

### Note on statistical power

With n = 3 per group, few or no proteins pass BH correction. The script
anticipates this: the volcano colours points by nominal p < 0.05 and
|log2FC| > 1, draws the BH FDR = 0.05 threshold as a dotted line, and if fewer
than 10 proteins clear FDR it labels the strongest nominal hits instead. The
heatmap likewise falls back to the top 50 proteins by p-value, and says so in
its title. Read the panels as the script labels them: nominal unless stated
otherwise.

## Running this section

    bash 03_protein_corona_proteomics/run_section.sh

or the script directly, from any working directory:

    python3 03_protein_corona_proteomics/ms_analysis.py

| | |
|---|---|
| Command | `python3 03_protein_corona_proteomics/ms_analysis.py` (no arguments, no environment variables) |
| Reads | `data/03_protein_corona_proteomics/Riccado_Report.csv` |
| Writes | `results/03_protein_corona_proteomics/TBI_blast_proteomics_results.xlsx`, `results/03_protein_corona_proteomics/summary.json`, `results/03_protein_corona_proteomics/figures/*.{png,pdf,svg}` |
| Log | `results/03_protein_corona_proteomics/logs/ms_analysis.log` (runner only) |
| Run time | about 5 minutes (281–369 s over three runs on the machine used for the paper), most of it waiting for the web services and the three GSEA runs |
| Hardware | ordinary CPU machine; the GSEA step uses 4 threads |

`NANOAI_DATA` and `NANOAI_RESULTS` move the input and output roots (default `<repo>/data` and
`<repo>/results`). `data/` is never written to; `results/` is not tracked by git.

The runner first checks the input file, the Python packages and that the three services answer
(`preflight`), then runs the script, and prints one `PASS` / `FAIL` line for each with a non-zero
exit status if either failed.

### Network access is required

`ms_analysis.py` cannot run offline. No annotation or gene-set file is shipped or read from
disk; it queries three public services while it runs:

| Service | Endpoint | Used for |
|---|---|---|
| UniProt REST API | `https://rest.uniprot.org/uniprotkb/search` | primary gene symbol of each representative accession (stage 1), 80 accessions per request |
| g:Profiler g:Convert | `https://biit.cs.ut.ee/gprofiler/api/convert/convert/` | Ensembl gene identifier of each accession (stage 1) |
| g:Profiler g:GOSt | `https://biit.cs.ut.ee/gprofiler/api/gost/profile/` | over-representation tests (stage 6): sources GO:BP, GO:CC, GO:MF, KEGG, REAC, WP, CORUM, HPA; g:SCS threshold 0.05 |
| Enrichr, through `gseapy` | `https://maayanlab.cloud/Enrichr/` | gene-set libraries of the pre-ranked GSEA (stage 7): `GO_Biological_Process_2025`, `Reactome_Pathways_2024`, `KEGG_2021_Human`. `gseapy` caches them in `~/.cache/gseapy/` after the first download |

Only protein accessions and gene identifiers from the report are sent.

Because the annotation comes from live databases, the gene symbols and the enrichment tables
follow whatever release the services serve on the day of the run. The differential-abundance
statistics, the PCA and the counts do not depend on them.

### Python packages

`numpy`, `pandas`, `scipy`, `statsmodels`, `matplotlib`, `requests`, `adjustText`, `gseapy` and
`openpyxl`. Versions used: Python 3.13, NumPy 2.5.1, pandas 3.0.3, SciPy 1.18.0,
statsmodels 0.15.0, matplotlib 3.11.0, adjustText 1.4.0, openpyxl 3.1.5, requests 2.34.2,
gseapy 1.3.1.

### Outputs

| File | Content |
|---|---|
| `summary.json` | counts: protein groups in the report, detected, tested, differentially abundant, nominal p < 0.05, presence/absence |
| `TBI_blast_proteomics_results.xlsx` | sheets `Summary`, `QC_runs`, `1_Detected_proteins`, `1_Enrich_detected_vs_genome`, `1_Enrich_top10pct_abundant`, `2_DA_blast_vs_ctrl`, `2_DAPs`, `2_On-Off_proteins`, `3_Enrich_DAP_*`, `3_Enrich_relaxed_*`, `3_GSEA_preranked` |
| `figures/01_QC_proteins_per_run`, `02_sample_correlation` | quality control |
| `figures/03_PCA` | PCA of the proteins quantified in all six samples |
| `figures/04_volcano` | volcano plot |
| `figures/05_DAP_heatmap` | top-protein heatmap |
| `figures/06_enrich_detected_proteome`, `07_enrich_top10pct_abundant` | enrichment dot plots |
| `figures/11_GSEA` | pre-ranked GSEA, terms with FDR q < 0.25 |

Figures `08`–`10c` (enrichment of the differentially abundant sets) are written only when
g:Profiler returns significant terms for those sets; with the published data it returns none.

### What a rerun reproduces

Checked on the machine used for the paper, same package versions, against the workbook of the
original run: `summary.json` is identical, and every numeric column of every sheet except
`3_GSEA_preranked` is identical (differential abundance, detected proteome, both enrichment
tables including adjusted p values and fold enrichments). What differed is what
the live databases had changed in between: UniProt returned a different primary gene symbol for
5 of the 6,246 representative accessions, so those symbols differ in the `Gene` column and in the
`genes` lists of the enrichment tables, and, because the GSEA ranks proteins by gene symbol, 19
of the 3,088 GSEA terms changed their enrichment score (ES and NES) and one term dropped below
the minimum set size. The FDR q value is computed across all terms of a library, so it shifted
for most terms (by at most 0.14). The nine GSEA terms with FDR q < 0.05 were the same. Two runs
on the same day gave identical workbooks and `summary.json`.

One figure is not identical from run to run: in `figures/04_volcano` the gene labels are placed
by `adjustText`, whose optimisation stops on a wall-clock limit, so some labels land in slightly
different positions each time (15 of 51 between two runs here). The plotted points and the set of
labelled genes are identical between runs and to the original figure.

## Inputs

| File | What it is |
|---|---|
| `../data/03_protein_corona_proteomics/Riccado_Report.csv` | Spectronaut 17 directDIA protein-group report. |

Requires `requests` (two g:Profiler REST endpoints), `statsmodels`, `scipy`,
`pandas`, `matplotlib` and `adjustText`. It also imports `gseapy` (stage 7) and writes the
workbook with `openpyxl`; see [Running this section](#running-this-section) for the full list
and for the network access the script needs.

Source data for this section is in [`../data/03_protein_corona_proteomics`](../data/03_protein_corona_proteomics); see
[`../data/README.md`](../data/README.md) for formats and column schemas.
