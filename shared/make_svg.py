"""
Re-render paper figures as editable SVG by running the ORIGINAL analysis scripts unchanged
and intercepting Figure.savefig. For every figure a script saves we write:
  <name>__full.svg          the whole figure
  <name>__panel_<x>.svg     each axes on its own (other axes/suptitle hidden, tight bbox)
Text stays as <text> (svg.fonttype='none') so it can be edited/re-sized in Illustrator,
Inkscape or PowerPoint; font is declared as Arial.

usage: python make_svg.py <script.py> <cwd> <out_subdir> [ENV=VAL ...]
"""
import os, sys, runpy, string
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

script, cwd, sub = sys.argv[1], sys.argv[2], sys.argv[3]
for kv in sys.argv[4:]:
    k, v = kv.split("=", 1); os.environ[k] = v
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), sub)
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

script = os.path.abspath(script)
sys.path.insert(0, os.path.dirname(script))
os.chdir(cwd)
sys.argv = [script]
runpy.run_path(script, run_name="__main__")
