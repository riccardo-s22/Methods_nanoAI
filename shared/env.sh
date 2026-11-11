#!/usr/bin/env bash
# Portable path resolution for the shell scripts. Source it, do not execute it:
#
#     source "$(dirname "${BASH_SOURCE[0]}")/../shared/env.sh"
#
# It exports the same three locations shared/repo_paths.py resolves, so shell
# and Python agree:
#
#     NANOAI_REPO      repository root
#     NANOAI_DATA      source data      (default <repo>/data)
#     NANOAI_RESULTS   outputs          (default <repo>/results)
#
# NANOAI_DATA and NANOAI_RESULTS are respected if already set.

NANOAI_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export NANOAI_REPO
export NANOAI_DATA="${NANOAI_DATA:-${NANOAI_REPO}/data}"
export NANOAI_RESULTS="${NANOAI_RESULTS:-${NANOAI_REPO}/results}"
export PYTHON="${PYTHON:-python3}"
