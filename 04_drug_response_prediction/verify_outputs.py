"""
Check a rerun of this section against the published values. Not an analysis: it only reads.

  1. Tables: every regenerated CSV under results/ is compared, cell by cell, with the reference copy
     shipped in data/04_drug_response_prediction/ (row order, labels, max absolute difference).
  2. Statistics: the lines of the run logs that carry the reported numbers are compared with the
     strings the published runs printed.

Each check prints OK, MISMATCH or NOT RUN (the step that produces it was skipped). Exit status is
non-zero if anything is a MISMATCH. Run it after run_section.sh; run_section.sh calls it last.

    python 04_drug_response_prediction/verify_outputs.py
"""
import os, re, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from repo_paths import data_path, RESULTS_DIR
import numpy as np, pandas as pd

SEC = "04_drug_response_prediction"
RES = os.path.join(RESULTS_DIR, SEC)
LOGS = os.path.join(RES, "logs")
# The tables are expected to be identical to the shipped ones (they are on the machine the paper
# was run on). The tolerances only absorb last-digit differences between numerical libraries:
# one unit in the fifth decimal of the score tables, 1e-6 relative in the feature table
# (whose values are of order 1e-14 to 1e-7, so an absolute tolerance would be meaningless there).
ABS_TOL, REL_TOL = 2e-5, 1e-6
bad = 0

def report(status, what, detail=""):
    global bad
    bad += status == "MISMATCH"
    print("%-9s %s%s" % (status, what, ("  [" + detail + "]") if detail else ""))

# ---------------------------------------------------------------- 1. tables
TABLES = [
    ("exp3_3H_organoid_scores.csv",
     os.path.join("per_well_scores", "exp3_3H_organoid_scores.csv"), "abs"),
    (os.path.join("sensitivity_true_log_ratios", "exp3_3H_organoid_scores.csv"),
     os.path.join("per_well_scores", "exp3_3H_organoid_scores_truelog.csv"), "abs"),
    ("exp3_3H_chirality_peak_features.csv",
     os.path.join("chirality_features", "exp3_3H_chirality_peak_features.csv"), "rel"),
]
print("== tables: results/%s/<file>  vs  data/%s/<file>" % (SEC, SEC))
for new, ref, mode in TABLES:
    pn, pr_ = os.path.join(RES, new), data_path(SEC, ref)
    if not os.path.exists(pn):
        report("NOT RUN", new); continue
    a, b = pd.read_csv(pn), pd.read_csv(pr_)
    if a.shape != b.shape or list(a.columns) != list(b.columns):
        report("MISMATCH", new, "shape %s vs %s, columns equal: %s" % (a.shape, b.shape, list(a.columns) == list(b.columns)))
        continue
    num = [c for c in a.columns if pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c])]
    lab = [c for c in a.columns if c not in num]
    labels_ok = all((a[c].astype(str) == b[c].astype(str)).all() for c in lab)
    A, B = a[num].to_numpy(float), b[num].to_numpy(float)
    d = np.abs(A - B)
    mx = float(np.nanmax(d)) if d.size else 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(d == 0, 0.0, d / np.abs(B))
    mxrel = float(np.nanmax(rel)) if d.size else 0.0
    ok = labels_ok and (mx <= ABS_TOL if mode == "abs" else mxrel <= REL_TOL)
    report("OK" if ok else "MISMATCH", new,
           "%d rows x %d cols, row order and labels identical: %s, max abs diff %.3g, max rel diff %.3g, identical cells %d/%d"
           % (a.shape[0], a.shape[1], labels_ok, mx, mxrel, int((A == B).sum()), A.size))

# ---------------------------------------------------------------- 2. statistics printed in the logs
# (log file, what it is, the text the published run printed; whitespace runs are collapsed before comparing)
CHECKS = [
    ("exp3_perm_power.log", "Supp. Fig. 7A, shape vectors: observed AUC", "observed LOO-LDA AUC = 0.9033 (PCA 20-dim, n=32/32)"),
    ("exp3_perm_power.log", "Supp. Fig. 7A: null distribution", "permutation null (B0=20000): mean=0.474 sd=0.099 max=0.850 95th=0.636 99th=0.698"),
    ("exp3_perm_power.log", "Supp. Fig. 7A: z", "normal-tail: z=4.32 p_normal=7.69e-06g"),
    ("exp3_perm_power.log", "Supp. Fig. 7B, 66 ratios: observed AUC", "observed LOO-LDA AUC = 0.9375 (PCA 15-dim, n=32/32)"),
    ("exp3_perm_power.log", "Supp. Fig. 7B: null distribution", "permutation null (B0=20000): mean=0.473 sd=0.102 max=0.833 95th=0.635 99th=0.696"),
    ("exp3_perm_power.log", "Supp. Fig. 7B: z", "normal-tail: z=4.54 p_normal=2.76e-06g"),
    ("exp3_perm_power.log", "Supp. Fig. 7A: never exceeded in 20,000 permutations",
     "Method A: FULL-SHAPE (2560-d -> PCA20) AUC=0.903 z=4.3 -> effect is strong; null max over 20000 perms=0.850 (obs NEVER exceeded)"),
    ("exp3_perm_power.log", "Supp. Fig. 7B: never exceeded in 20,000 permutations",
     "Method B: FITTED-PEAK RATIOS (66 -> PCA15) AUC=0.938 z=4.5 -> effect is strong; null max over 20000 perms=0.833 (obs NEVER exceeded)"),
    ("exp3_drug_axis.log", "Fig. 4b: sham", "sham -0.220 +/- 0.111"),
    ("exp3_drug_axis.log", "Fig. 4b: blast", "tbi 0.220 +/- 0.124"),
    ("exp3_drug_axis.log", "Fig. 4b: blast + resveratrol", "resv_tbi 0.029 +/- 0.187"),
    ("exp3_drug_axis.log", "Fig. 4b: blast + MG132", "mg132_tbi 0.199 +/- 0.133"),
    ("exp3_drug_axis.log", "Fig. 4b: sham + resveratrol", "resv_sham -0.370 +/- 0.207"),
    ("exp3_drug_axis.log", "Fig. 4b: sham + MG132", "mg132_sham -0.093 +/- 0.092"),
    ("exp3_drug_axis.log", "blast plate, 3-group permutation F", "TBI plate 3-group on axis: F=14.93 p=0.0001"),
    ("exp3_drug_axis.log", "blast vs blast + resveratrol", "tbi vs resv_tbi axis-diff=-0.191 p=0.0001"),
    ("exp3_drug_axis.log", "blast vs blast + MG132", "tbi vs mg132_tbi axis-diff=-0.021 p=0.5267"),
    ("exp3_drug_axis.log", "sham plate, 3-group permutation F", "SHAM plate 3-group on axis: F=23.03 p=0.0001"),
    ("exp3_drug_axis.log", "DiD resveratrol", "resveratrol DiD=-0.041 (SE 0.058) p=0.4796"),
    ("exp3_drug_axis.log", "DiD MG132", "MG132 DiD=-0.148 (SE 0.043) p=0.0008"),
    ("exp3_drug_axis.log", "resveratrol rescue", "resveratrol % rescue ALONG AXIS = 43% (95% CI 28..60%)"),
    ("exp3_drug_axis.log", "MG132 move", "MG132 % move ALONG AXIS = 5% (95% CI -10..18%)"),
    ("exp3_drug_axis.log", "Supp. Fig. 8B: resveratrol on blast", "resv on TBI 4.04 -0.19 4.03 -0.05"),
    ("exp3_drug_axis.log", "Supp. Fig. 8B: MG132 on blast", "MG132 on TBI 1.29 -0.02 1.29 -0.02"),
    ("exp3_drug_axis.log", "Supp. Fig. 8B: resveratrol on sham", "resv on sham 4.34 -0.15 4.33 -0.03"),
    ("exp3_drug_axis.log", "Supp. Fig. 8B: MG132 on sham", "MG132 on sham 3.84 +0.13 3.84 +0.03"),
    ("exp3_3H_organoid_score.log", "per-well score, LOO AUC", "LOO AUC sham vs TBI=0.949"),
    ("plot_wb_swcnt_corr.log", "Supp. Fig. 10b: NfL", "NfL: Pearson r=-0.809 p=0.0514 | Spearman rho=-0.829 p=0.0416"),
    ("plot_wb_swcnt_corr.log", "Supp. Fig. 10c: NfH", "NfH: Pearson r=-0.513 p=0.298 | Spearman rho=-0.543 p=0.266"),
    ("exp3_3H_pressure.log", "pressure test: shape AUC", "shape AUC: headline=0.903 | outliers-removed=0.907 | detrended=0.978 | detrended+clean=0.969"),
    # sensitivity re-run with true log ratios
    ("exp3_perm_power_truelog.log", "true log, 66 ratios: observed AUC", "observed LOO-LDA AUC = 0.8867 (PCA 15-dim, n=32/32)"),
    ("exp3_perm_power_truelog.log", "true log, 66 ratios: z", "normal-tail: z=4.02 p_normal=2.93e-05g"),
    ("exp3_perm_power_truelog.log", "true log, 66 ratios: never exceeded in 20,000 permutations",
     "Method B: FITTED-PEAK RATIOS (66 -> PCA15) AUC=0.887 z=4.0 -> effect is strong; null max over 20000 perms=0.811 (obs NEVER exceeded)"),
    ("exp3_perm_power_truelog.log", "true log: shape vectors unaffected", "observed LOO-LDA AUC = 0.9033 (PCA 20-dim, n=32/32)"),
    ("exp3_3H_organoid_score_truelog.log", "true log: per-well score, LOO AUC", "LOO AUC sham vs TBI=0.913"),
    ("exp3_drug_axis_truelog.log", "true log: blast vs blast + resveratrol", "tbi vs resv_tbi axis-diff=-0.711 p=0.0011"),
    ("exp3_drug_axis_truelog.log", "true log: resveratrol rescue", "resveratrol % rescue ALONG AXIS = 31% (95% CI 15..46%)"),
    ("exp3_drug_axis_truelog.log", "true log: DiD resveratrol", "resveratrol DiD=-0.078 (SE 0.287) p=0.7860"),
    ("exp3_drug_axis_truelog.log", "true log: DiD MG132", "MG132 DiD=-0.723 (SE 0.314) p=0.0229"),
    ("exp3_peakfit_ratios_truelog.log", "true log: 12 fractions", "[12 AUC-fractions ] sham vs TBI LOO-LDA AUC=0.613 perm-p=0.0765"),
    ("exp3_peakfit_ratios_truelog.log", "true log: 66 ratios", "[66 log-ratios ] sham vs TBI LOO-LDA AUC=0.887 perm-p=0.0005"),
    ("exp3_peakfit_ratios_truelog.log", "true log: 78 combined", "[both (78) ] sham vs TBI LOO-LDA AUC=0.895 perm-p=0.0005"),
    ("exp3_peakfit_ratios_truelog.log", "true log: ratios with q<0.05", "# ratios with q<0.05: 17 / 66"),
    ("exp3_peakfit_ratios_truelog__FLO1e-6.log", "FLO=1e-6 control: 66 ratios (published)", "[66 log-ratios ] sham vs TBI LOO-LDA AUC=0.938 perm-p=0.0005"),
    ("exp3_peakfit_ratios_truelog__FLO1e-6.log", "FLO=1e-6 control: 78 combined (published)", "[both (78) ] sham vs TBI LOO-LDA AUC=0.920 perm-p=0.0005"),
    ("exp3_peakfit_ratios_truelog__FLO1e-6.log", "FLO=1e-6 control: ratios with q<0.05", "# ratios with q<0.05: 0 / 66"),
    ("exp3_peakfit_ratios_truelog__FLO1e-6.log", "FLO=1e-6 control: top ratio", "(8,3)/(10,2) d=-0.86 AUC=0.75 p=7.7e-04 q=0.05"),
    ("plot_wb_swcnt_corr_truelog.log", "true log: NfL", "NfL: Pearson r=-0.753 p=0.084 | Spearman rho=-0.771 p=0.0724"),
    ("plot_wb_swcnt_corr_truelog.log", "true log: NfH", "NfH: Pearson r=-0.438 p=0.385 | Spearman rho=-0.371 p=0.468"),
]
squash = lambda s: re.sub(r"\s+", " ", s)
print("\n== statistics: results/%s/logs/<log>" % SEC)
cache = {}
for log, what, expect in CHECKS:
    p = os.path.join(LOGS, log)
    if not os.path.exists(p):
        report("NOT RUN", "%s: %s" % (log, what)); continue
    if log not in cache:
        cache[log] = squash(open(p, encoding="utf-8", errors="replace").read())
    report("OK" if squash(expect) in cache[log] else "MISMATCH", "%s: %s" % (log, what), expect)

print("\n%s" % ("all checks that ran agree with the published values" if not bad else "%d MISMATCH" % bad))
sys.exit(1 if bad else 0)
