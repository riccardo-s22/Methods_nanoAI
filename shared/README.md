# Shared tooling

`make_svg.py` regenerates the paper figures as true vector SVG by re-running the
original analysis scripts unchanged, with a Matplotlib SVG backend, and splitting
each figure into its panels.

    python make_svg.py <script.py> <working dir> <out subdir> [ENV=VAL ...]

Outputs are named `*__full.svg` (the whole figure) and `*__panel_X.svg` (each
subplot alone, lettered in reading order, left to right then top to bottom).
Text stays editable; the font is Arial.

Per-section invocation:

    section 1   cwd = exp2/   EXP2_DIR=exp2/{2h,24h} EXP2_TAG={2h,24h} MIX_INVERTED=1
    section 4   cwd = the directory holding exp3_3H_peakfeats.npz   EXP3_TP=3H
