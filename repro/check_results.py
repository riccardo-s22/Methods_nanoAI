"""
Compare the tables a rerun wrote under results/ with the reference tables shipped in data/.

    python3 repro/check_results.py

It reads data files only (CSV / TSV), never log text, and changes nothing. bash reproduce.sh runs
it last. One line per table:

    PASS          same shape, same columns, same labels, numbers equal within the tolerance
    FAIL          anything else; the line says what differs
    NOT PRODUCED  the regenerated table does not exist (its section or step was not run)

Tolerances follow the precision the reference file is printed with:

    per-well score tables       5 decimals printed       |a - b| <= 2e-5
    every other table           full double precision    |a - b| <= 1e-6 * |b|   (values go down to 1e-14)

Labels must be equal as text, except that characters outside ASCII are ignored (the shipped FT
pathway matrix has "IKK?" where Enrichr returns the Greek letter).

What a PASS means differs between rows; the "source" printed with each row says which:

    recomputed    the table is computed again from the spectra / the shipped per-gene tables
    pass-through  in a default run section 02 reads the two Supplementary Figure 6 matrices from
                  data/ and writes them back, so the row only shows the round trip is lossless.
                  They are recomputed (Enrichr) only by run_section.sh --with-project-tree.

Tables without a reference in data/ are not checked here: docs/EXPECTED_RESULTS.md lists where to
read their numbers. Exit status: 1 if any row is FAIL, otherwise 0.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from repo_paths import DATA_DIR, RESULTS_DIR  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

S2 = "02_geneformer_in_silico_screen"
S4 = "04_drug_response_prediction"

# (regenerated file under results/, reference file under data/, tolerance mode, source, produced by)
TABLES = [
    (f"{S4}/exp3_3H_organoid_scores.csv",
     f"{S4}/per_well_scores/exp3_3H_organoid_scores.csv", "abs", "recomputed",
     "exp3_3H_organoid_score.py"),
    (f"{S4}/sensitivity_true_log_ratios/exp3_3H_organoid_scores.csv",
     f"{S4}/per_well_scores/exp3_3H_organoid_scores_truelog.csv", "abs", "recomputed",
     "exp3_3H_organoid_score_truelog.py"),
    (f"{S4}/exp3_3H_chirality_peak_features.csv",
     f"{S4}/chirality_features/exp3_3H_chirality_peak_features.csv", "rel", "recomputed",
     "exp3_export_peakfeats.py"),
    (f"{S2}/drug_repurpose/crossmodel_consensus.tsv",
     f"{S2}/drug_signature_enrichment/crossmodel_consensus.tsv", "rel", "recomputed",
     "excit_drug_crossmodel.py"),
    (f"{S2}/pathway_heatmap/excit_pathway_enrichment_matrix_FT.tsv",
     f"{S2}/pathway_enrichment/excit_pathway_enrichment_matrix_FT.tsv", "rel", "pass-through in a default run",
     "excit_pathway_heatmap_FT.py"),
    (f"{S2}/pathway_heatmap/excit_pathway_enrichment_matrix_MT.tsv",
     f"{S2}/pathway_enrichment/excit_pathway_enrichment_matrix_MT.tsv", "rel", "pass-through in a default run",
     "excit_pathway_heatmap_MT.py"),
    # written only by run_section.sh --with-project-tree
    (f"{S2}/drug_repurpose/FT/drug_candidates_ranked_union.tsv",
     f"{S2}/drug_signature_enrichment/FT/drug_candidates_ranked_union.tsv", "rel", "recomputed",
     "excit_drug_repurpose.py FT union (--with-project-tree)"),
    (f"{S2}/drug_repurpose/FT/drug_candidates_ranked_core.tsv",
     f"{S2}/drug_signature_enrichment/FT/drug_candidates_ranked_core.tsv", "rel", "recomputed",
     "excit_drug_repurpose.py FT core (--with-project-tree)"),
    (f"{S2}/drug_repurpose/FT/all_sig_drug_terms.tsv",
     f"{S2}/drug_signature_enrichment/FT/all_sig_drug_terms.tsv", "rel", "recomputed",
     "excit_drug_repurpose.py FT core (--with-project-tree)"),
]
ABS_TOL = 2e-5
REL_TOL = 1e-6


def read(path):
    return pd.read_csv(path, sep="\t" if path.endswith(".tsv") else ",")


def ascii_only(series):
    """Labels as text, with missing values and non-ASCII characters made comparable."""
    text = series.astype(object).where(series.notna(), "").astype(str)
    return text.str.replace(r"[^\x00-\x7f]", "?", regex=True)


def compare(new_path, ref_path, mode):
    """Returns (ok, detail)."""
    new, ref = read(new_path), read(ref_path)
    if list(new.columns) != list(ref.columns):
        return False, "columns differ: %s vs %s" % (list(new.columns)[:8], list(ref.columns)[:8])
    if new.shape != ref.shape:
        return False, "shape %s, reference %s" % (new.shape, ref.shape)
    numeric = [c for c in ref.columns
               if pd.api.types.is_numeric_dtype(ref[c]) and pd.api.types.is_numeric_dtype(new[c])
               and not pd.api.types.is_bool_dtype(ref[c])]
    labels = [c for c in ref.columns if c not in numeric]
    for c in labels:
        differ = int((ascii_only(new[c]) != ascii_only(ref[c])).sum())
        if differ:
            return False, "%d of %d labels differ in column %r" % (differ, len(ref), c)
    worst_abs = worst_rel = 0.0
    if numeric:
        a = new[numeric].to_numpy(float)
        b = ref[numeric].to_numpy(float)
        if not np.array_equal(np.isnan(a), np.isnan(b)):
            return False, "missing values are not in the same cells"
        if not np.array_equal(np.isinf(a), np.isinf(b)):
            return False, "infinite values are not in the same cells"
        finite = np.isfinite(a) & np.isfinite(b)
        if not np.array_equal(a[~finite & ~np.isnan(a)], b[~finite & ~np.isnan(b)]):
            return False, "infinite values differ in sign"
        d = np.abs(a[finite] - b[finite])
        scale = np.abs(b[finite])
        if d.size:
            worst_abs = float(d.max())
            with np.errstate(divide="ignore", invalid="ignore"):
                rel = np.where(d == 0, 0.0, d / scale)
            worst_rel = float(rel.max())
            limit = ABS_TOL if mode == "abs" else REL_TOL * scale
            n_out = int((d > limit).sum())
            if n_out:
                return False, "%d of %d numbers outside tolerance (max abs diff %.3g, max rel diff %.3g)" % (
                    n_out, d.size, worst_abs, worst_rel)
    tol = "abs tol %g" % ABS_TOL if mode == "abs" else "rel tol %g" % REL_TOL
    return True, "%d x %d, %d label and %d numeric columns, max abs diff %.3g, %s" % (
        ref.shape[0], ref.shape[1], len(labels), len(numeric), worst_abs, tol)


def main():
    print("regenerated tables in  %s" % RESULTS_DIR)
    print("reference tables in    %s" % DATA_DIR)
    counts = {"PASS": 0, "FAIL": 0, "NOT PRODUCED": 0}
    for new, ref, mode, source, script in TABLES:
        new_path = os.path.join(RESULTS_DIR, *new.split("/"))
        ref_path = os.path.join(DATA_DIR, *ref.split("/"))
        if not os.path.isfile(ref_path):
            status, detail = "FAIL", "reference missing: data/%s" % ref
        elif not os.path.isfile(new_path):
            status, detail = "NOT PRODUCED", "written by %s" % script
        else:
            try:
                ok, detail = compare(new_path, ref_path, mode)
            except Exception as ex:  # unreadable or malformed table
                ok, detail = False, "could not compare: %s" % ex
            status = "PASS" if ok else "FAIL"
            detail = "%s; %s" % (source, detail)
        counts[status] += 1
        print("%-12s  results/%s\n              vs data/%s\n              %s" % (status, new, ref, detail))
    print("check_results: %d PASS, %d FAIL, %d NOT PRODUCED" % (
        counts["PASS"], counts["FAIL"], counts["NOT PRODUCED"]))
    return 1 if counts["FAIL"] else 0


if __name__ == "__main__":
    sys.exit(main())
