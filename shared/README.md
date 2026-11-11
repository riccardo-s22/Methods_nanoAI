# Shared tooling

## Path resolution

`repo_paths.py` (Python) and `env.sh` (shell) resolve every input and output location of the
repository, so no script carries an absolute path and all of them run from any working directory:

    inputs    <repo>/data/...       override with NANOAI_DATA
    outputs   <repo>/results/...    override with NANOAI_RESULTS

`data/` is read-only; everything a rerun writes goes under `results/`, which git ignores.

`repo_paths.py` holds locations only, never an analysis parameter. Section 02 adds its own
project roots on top of it (`02_geneformer_in_silico_screen/gf_paths.py`, `gf_env.sh`).
The top-level `reproduce.sh` and every `<section>/run_section.sh` source `env.sh`, so setting
`NANOAI_DATA`, `NANOAI_RESULTS` or `PYTHON` in the shell applies to a whole rerun.

## make_svg.py

`make_svg.py` regenerates the paper figures as true vector SVG by re-running the
original analysis scripts unchanged, with a Matplotlib SVG backend, and splitting
each figure into its panels.

    python3 shared/make_svg.py <script.py> [<working dir>] [<out subdir>] [ENV=VAL ...]

- `<script.py>`: absolute, relative to the current directory, or relative to the repository
  root (for example `01_swcnt_blast_classification/exp2_analysis2.py`).
- `<out subdir>`: where the SVGs go, relative to the results tree (`<repo>/results`, or
  `$NANOAI_RESULTS`); an absolute path is used as given. Default: `<section of the script>/svg`.
- `<working dir>`: optional and only accepted together with `<out subdir>`; the scripts no
  longer depend on the working directory, so it can be left out.
- `ENV=VAL`: environment variables set for the script.

Outputs are named `*__full.svg` (the whole figure) and `*__panel_X.svg` (each
subplot alone, lettered in reading order, left to right then top to bottom).
Text stays editable; the font is Arial.

The harness intercepts `savefig`, so during this run the script's own PNG or PDF files are not
written; anything else the script writes (tables, JSON, text) goes where it normally does. Each
call runs the whole script again, including its statistics.

Set `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` in the shell for these calls.
The figure scripts set them at their top, but the harness has already imported NumPy by then, so
without the shell setting BLAS runs multi-threaded and the permutation-heavy scripts take many
times longer. Passing them as `ENV=VAL` arguments is also too late.

Per-section invocation, from the repository root:

    section 1   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 shared/make_svg.py 01_swcnt_blast_classification/exp2_analysis2.py 01_swcnt_blast_classification/svg/exp2_24h EXP2_TAG=24h EXP2_DIR=data/01_swcnt_blast_classification/spectra/24h MIX_INVERTED=1
                (EXP2_TAG and the EXP2_DIR folder are 2h or 24h, set together; done for the published panels by 01_swcnt_blast_classification/run_section.sh)
    section 4   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 shared/make_svg.py 04_drug_response_prediction/exp3_drug_axis.py 04_drug_response_prediction/svg EXP3_TP=3H
                (done for the published panels by 04_drug_response_prediction/run_section.sh, which also sets the thread variables)
