"""Portable path resolution for every script in this repository.

The scripts were first run with absolute paths of the machines they ran on.
They now resolve inputs and outputs through this module instead, so a fresh
clone runs without editing anything:

    inputs   <repo>/data/...       override with NANOAI_DATA
    outputs  <repo>/results/...    override with NANOAI_RESULTS

Only path resolution lives here. No analysis parameter is read from this file.

Usage, from a script N directories below the repository root:

    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    *[".."] * N, "shared"))
    from repo_paths import data_path, results_dir

    spectra = data_path("01_swcnt_blast_classification", "spectra", "24h")
    out = results_dir("01_swcnt_blast_classification")

Both helpers return plain strings, because the scripts build paths with
os.path.join and string formatting.
"""
import os

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.abspath(os.environ.get("NANOAI_DATA") or os.path.join(REPO_ROOT, "data"))
RESULTS_DIR = os.path.abspath(os.environ.get("NANOAI_RESULTS") or os.path.join(REPO_ROOT, "results"))


def repo_path(*parts):
    """Path inside the repository. Not created, not checked."""
    return os.path.join(REPO_ROOT, *parts)


def data_path(*parts):
    """Path inside the read-only source-data tree. Not created, not checked."""
    return os.path.join(DATA_DIR, *parts)


def results_dir(*parts):
    """Output directory under the results tree, created if missing."""
    path = os.path.join(RESULTS_DIR, *parts)
    os.makedirs(path, exist_ok=True)
    return path


def env_path(name, default):
    """An existing per-script override (EXP2_DIR, EXP2_OUT, EXP3_DIR, ...) if
    set and non-empty, otherwise the repository default."""
    return os.environ.get(name) or default
