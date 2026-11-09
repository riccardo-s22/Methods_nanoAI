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
(`../Riccado_Report.csv`); outputs are the figure panels plus
`TBI_blast_proteomics_results.xlsx` and `summary.json`.

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

## Inputs

| File | What it is |
|---|---|
| `Riccado_Report.csv` | Spectronaut 17 directDIA protein-group report. Not in this repository; available from the corresponding author. |

Requires `requests` (two g:Profiler REST endpoints), `statsmodels`, `scipy`,
`pandas`, `matplotlib` and `adjustText`.

Source data for this section is in [`../data/03_protein_corona_proteomics`](../data/03_protein_corona_proteomics); see
[`../data/README.md`](../data/README.md) for formats and column schemas.
