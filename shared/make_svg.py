"""
Re-render paper figures as editable SVG by running the ORIGINAL analysis scripts unchanged
and intercepting Figure.savefig. For every figure a script saves we write:
  <name>__full.svg          the whole figure
  <name>__panel_<x>.svg     each axes on its own (other axes/suptitle hidden, tight bbox)
Text stays as <text> (svg.fonttype='none') so it can be edited/re-sized in Illustrator,
Inkscape or PowerPoint; font is declared as Arial.

usage: python make_svg.py <script.py> [<cwd>] [<out_subdir>] [ENV=VAL ...]

  <script.py>   the figure script: an absolute path, a path relative to the current directory,
                or a path relative to the repository root
                (e.g. 01_swcnt_blast_classification/exp2_analysis2.py)
  <cwd>         optional, only together with <out_subdir>: directory to run the script in.
                The scripts no longer depend on the working directory, so it can be left out.
  <out_subdir>  optional: where the SVGs go, relative to the results tree
                (<repo>/results, or $NANOAI_RESULTS); an absolute path is used as given.
                Default: <section directory of the script>/svg
  ENV=VAL       environment variables for the script (EXP2_TAG=24h, MIX_INVERTED=1, EXP3_TP=3H, ...)

Because savefig is intercepted, the script's own PNG/PDF files are not written by this run;
anything else it writes (tables, JSON, text) goes where the script normally puts it.
"""
import os, sys, runpy, string
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

args = sys.argv[1:]
npos = next((i for i, a in enumerate(args) if "=" in a), len(args))   # positionals come before ENV=VAL
pos, kvs = args[:npos], args[npos:]
if not 1 <= len(pos) <= 3:
    sys.exit(__doc__)
script = pos[0]
cwd = pos[1] if len(pos) == 3 else None
sub = pos[-1] if len(pos) >= 2 else None
for kv in kvs:
    k, v = kv.split("=", 1); os.environ[k] = v

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))        # after ENV=VAL, so NANOAI_* given there count
from repo_paths import REPO_ROOT, RESULTS_DIR

if not os.path.exists(script) and os.path.exists(os.path.join(REPO_ROOT, script)):
    script = os.path.join(REPO_ROOT, script)                          # given relative to the repository root
script = os.path.abspath(script)
if sub is None:                                                       # default: <section>/svg under the results tree
    rel = os.path.relpath(script, REPO_ROOT)
    sub = os.path.join(rel.split(os.sep)[0], "svg") if os.sep in rel and not rel.startswith("..") else "svg"
OUTDIR = sub if os.path.isabs(sub) else os.path.join(RESULTS_DIR, sub)
os.makedirs(OUTDIR, exist_ok=True)

plt.rcParams.update({"svg.fonttype": "none", "font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "DejaVu Sans"]})

def _svg_save(fig, fname, *a, **kw):
    base = os.path.splitext(os.path.basename(str(fname)))[0]
    fig.savefig_orig(os.path.join(OUTDIR, base + "__full.svg"), format="svg", bbox_inches="tight")
    axes = [ax for ax in fig.axes if ax.get_visible()]
    texts = list(fig.texts) + ([fig._suptitle] if fig._suptitle else [])
    legends = list(fig.legends)
    vis = {id(o): o.get_visible() for o in axes + texts + legends}
    if len(axes) > 1:
        for t in texts + legends: t.set_visible(False)
        for i, ax in enumerate(axes):
            for o in axes: o.set_visible(o is ax)
            fig.savefig_orig(os.path.join(OUTDIR, "%s__panel_%s.svg" % (base, string.ascii_lowercase[i])),
                             format="svg", bbox_inches="tight")
        for o in axes + texts + legends: o.set_visible(vis[id(o)])
    print("[svg] %s -> %d panel(s)" % (base, len(axes)), flush=True)

Figure.savefig_orig = Figure.savefig
Figure.savefig = _svg_save

sys.path.insert(0, os.path.dirname(script))
if cwd is not None:
    os.chdir(cwd)
sys.argv = [script]
runpy.run_path(script, run_name="__main__")
