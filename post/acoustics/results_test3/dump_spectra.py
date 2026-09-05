#!/usr/bin/env python3
"""Robust ensemble line spectra for Gate A: per-k2 MEDIAN (and mean) of the
component spectral densities over all realizations, at each target strain.

Heavy-tailed intermittency makes the arithmetic ensemble mean non-convergent
even at 1024 realizations (see diag_outliers.py), so the fit targets are the
median spectra; the mean is kept alongside for comparison.  All realizations
share the FFT k2 grid at fixed strain (deterministic dilatation).

    python3 dump_spectra.py <root> [<root2> ...]

Writes spectra_<case>.npz next to this script: k2_e<j>, med_E1_e<j>,
med_E2_e<j>, med_E3_e<j>, mean_E*_e<j> for each strain index j, plus strains,
nrlz, L_e<j>.
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import odt_io                                          # noqa: E402

STRAINS = [0.0, 1.0, 2.0, 3.0, 3.9]
NU = 2048

# Override for cases whose dumps sit at other times (e.g. gateA_S1: precursor
# to t=0.4, dumps at t = e + 0.4):  DUMP_TIMES="0.4,0.9,1.4,1.9,2.4".
# load_fit_target matches on t (smag=1), so pass TIMES; "strains" in the npz
# then records these times — subtract the precursor offset downstream.
if os.environ.get("DUMP_TIMES"):
    STRAINS = [float(x) for x in os.environ["DUMP_TIMES"].split(",")]
if os.environ.get("DUMP_NU"):
    NU = int(os.environ["DUMP_NU"])

for root in sys.argv[1:]:
    dirs = sorted(glob.glob(os.path.join(root, "data", "data_*")))
    n = len(dirs)
    out = {"strains": np.array(STRAINS), "nrlz": n}
    for j, e in enumerate(STRAINS):
        E = None
        for i, d in enumerate(dirs):
            t = odt_io.load_fit_target(d, strain=e, Nu=NU)
            if E is None:
                out[f"k2_e{j}"] = t["k2"]
                out[f"L_e{j}"] = t["L"]
                E = np.empty((3, n, t["k2"].size))
            E[0, i] = t["E1"]
            E[1, i] = t["E2"]
            E[2, i] = t["E3"]
        for c in range(3):
            out[f"med_E{c+1}_e{j}"] = np.median(E[c], axis=0)
            out[f"mean_E{c+1}_e{j}"] = np.mean(E[c], axis=0)
    fn = os.path.join(os.path.dirname(__file__),
                      "spectra_" + os.path.basename(root.rstrip("/")) + ".npz")
    np.savez_compressed(fn, **out)
    print("wrote", fn, f"({n} rlz)")
