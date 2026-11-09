Supplementary figure: blast vs sham classification (exp2, 2 h and 24 h)
=======================================================================
EXPERIMENTAL UNIT = one organoid = one well.  n = 32 sham + 32 blast per time point.
Sham = plate columns 2,4,6,8; blast = columns 3,5,7,9 (interleaved, so any smooth
plate gradient is balanced between groups by design).

FILES
  supp_auc_blast_vs_sham.{svg,pdf,png}  the figure (a,b ROC; c held-out scores;
                                        d permutation null; e AUC by held-out unit)
  make_supp_auc_compute.py              step 1 -> supp_auc_results.json
  make_supp_auc_figure.py               step 2 -> the figure
  supp_auc_results.json                 ROC coords, AUCs, scores, 1000-perm nulls
  check_simple_descriptors.py           can t-tests replace the classifier? (see below)
  check_permutation_power.py            power analysis of the permutation tests

Pipelines are the ORIGINAL ones in TBI/exp2 (exp2_load.py; the PCA20 classifier
definitions of exp2_classifier_compare.py; the row-out scheme of exp2_analysis2.py),
with THREE corrections:

 1. PCA(20) -> PCA(20, random_state=0).
    On 64 x 2560 sklearn silently picks the RANDOMIZED svd solver; with
    random_state=None the subspace differs between runs.
 2. AUC from decision_function, not predict_proba.
    LDA probabilities saturate at exactly 1.0, so wells tie at the float ceiling and
    roc_auc_score charges those pairs 0.5.
    -> (1)+(2) made the 2 h row-out AUC print as 1.000 or 0.9844 at random (3/10 runs
       gave 0.9844). Both are still present in exp2_analysis2.py:25,32 and
       exp2_classifier_compare.py:56.
 3. Scaler+PCA refit INSIDE every CV fold for the row-out and column-out schemes
    (exp2_analysis2.py:25 fits PCA on all 64 wells before the CV loop).

RESULTS (all corrections applied; every value reproducible)
  AUC = 1.0000 for all 5 classifiers x 3 held-out units x 2 time points (30/30).
  LDA score gap (min blast - max sham), leave-one-organoid-out: +250 (2 h), +332 (24 h).
  Label permutation, 1000 shuffles, leave-one-organoid-out LDA:
     2 h   observed 1.000; null mean 0.467, 95th 0.636, max 0.747; p = 0.001
     24 h  observed 1.000; null mean 0.471, 95th 0.633, max 0.729; p = 0.001

WHY THE CLASSIFIER IS NEEDED (check_simple_descriptors.py)
  24 h: simple descriptors work. Emission centroid red-shifts +6.0 to +8.5 nm at all
        five excitations; best single descriptor AUC 0.858 (ex640, d=1.62, MWU
        p=8.6e-7); dilution series on mean centroid rho=+0.486, ANOVA p=5.8e-6.
   2 h: simple descriptors FAIL. Best centroid AUC 0.659 (ex750, p=0.029 uncorrected,
        0.145 after correcting for 5 excitations); mean centroid AUC 0.549 (p=0.51);
        total intensity AUC 0.588 (p=0.23); dilution trend on centroid rho=+0.066
        (p=0.54). The 2 h signal is multivariate spectral SHAPE with no univariate
        correlate, so the classifier is load-bearing at 2 h, not decorative.
