# vendored/

Modules the scripts of this section import that lived outside the repository
when the analyses were run. They are copied here so that a fresh clone runs.
Nothing else from their source folders is included.

| File | Original location | Imported by | What is used |
|---|---|---|---|
| `chirality_interpolated_extraction.py` | `physical_descriptors/` folder of the laboratory's SWCNT descriptor project (`exp5/physical_descriptors/`) | `exp3_peakfit.py`, `exp3_peakfit_ratios.py`, `exp3_build_cache.py`, `supporting_checks/exp3_3H_pressure.py`, `sensitivity_true_log_ratios/exp3_peakfit_ratios_truelog.py` | `CHIRALITY_POSITIONS`, `INTERP_FACTOR`, `SEARCH_RADIUS_NM`, `build_interpolated_surface`, `find_local_max_in_radius`, `extract_consistent_gradient_points`, `fit_gaussian_and_extract` |
| `exp3_7d_load.py` | the drug-experiment working folder (`exp3/`), next to `exp3_load.py` | `supporting_checks/exp3_3H_pressure.py` | `outlier_keep()` only |

`chirality_interpolated_extraction.py` is a byte-for-byte copy
(MD5 `72775a03d51bdb65b62e3ca42741d06b`). It is the peak localization, windowing
and split-Gaussian fitting module described in the Methods. It is loaded by file
path (`importlib`), not through `sys.path`.

`exp3_7d_load.py` is the loader of the 7-day read, a separate experiment that is
not reported in the paper and whose spectra are not in the repository. Only its
`outlier_keep()` function is used here. It is identical to the original
(MD5 of the original `fe44ab3ca9c30cce97056900cfc7a0f0`) except for the default
of its data folder, which was an absolute path of the workstation and now
resolves through `shared/repo_paths.py`:

    -import numpy as np, glob, os
    -
    -BASE = os.environ.get("EXP3_DIR", "<absolute path of the original exp3 folder>")
    +import numpy as np, glob, os, sys
    +sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "shared"))
    +from repo_paths import data_path, env_path
    +
    +BASE = env_path("EXP3_DIR", data_path("04_drug_response_prediction", "spectra"))

`BASE` is read only by its `load_plate()`, which this section never calls.

Under NumPy 2.5 the fit-failure fallback of
`chirality_interpolated_extraction.py` calls `numpy.trapz`, which no longer
exists; see *Known limitations* in the section README. The module is left as it
was.
