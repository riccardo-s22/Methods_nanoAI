#!/bin/bash
# Re-run of the chirality-ratio analyses with true log10 ratios (flo = 0). Originals untouched.
T=/mnt/c/Users/riccardo-s/Documents/CNT/TBI; R=$T/paper/svg/truelog_rerun; PY=/home/riccardo-s/miniforge3/bin/python
E5=/mnt/c/Users/riccardo-s/Documents/CNT/targetALS/exp5
export PYTHONPATH=$T/exp3 EXP3_TP=3H MIX_INVERTED=1
cd $R && $PY exp2_chirality_classify_truelog.py > log_exp2_chirality.txt 2>&1 &
for s in exp3_drug_axis exp3_perm_power exp3_3H_organoid_score exp3_peakfit_ratios; do
  (cd $E5 && $PY $R/${s}_truelog.py > $R/log_$s.txt 2>&1 && \
   cd $T/paper/svg && $PY make_svg.py $R/${s}_truelog.py $E5 truelog_rerun/svg EXP3_TP=3H > $R/log_svg_$s.txt 2>&1) &
done
wait; echo ALL DONE
