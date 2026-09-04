#!/usr/bin/env python3
"""Dump per-realization band means for all strains to bands_<case>.npz.

    python3 dump_bands.py <root> [<root2> ...]

Arrays: E2lo, Eplo, E2hi, Ephi, u2frac  (n_rlz x n_strain), names (n_rlz),
strains.  Written next to this script.
"""
import glob
import os
import sys

import numpy as np

trapezoid = getattr(np, "trapezoid", np.trapz)

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import odt_io                                          # noqa: E402

STRAINS = [0.0, 1.0, 2.0, 3.0, 3.9]
BAND_LO = (30.0, 100.0)
BAND_HI = (300.0, 800.0)
NU = 2048

for root in sys.argv[1:]:
    dirs = sorted(glob.glob(os.path.join(root, "data", "data_*")))
    n = len(dirs)
    out = {k: np.zeros((n, len(STRAINS))) for k in
           ("E2lo", "Eplo", "E2hi", "Ephi", "u2frac")}
    for i, d in enumerate(dirs):
        for j, e in enumerate(STRAINS):
            t = odt_io.load_fit_target(d, strain=e, Nu=NU)
            k2 = t["k2"]
            Ep = 0.5 * (t["E1"] + t["E3"])
            E2 = t["E2"]
            mlo = (k2 >= BAND_LO[0]) & (k2 <= BAND_LO[1])
            mhi = (k2 >= BAND_HI[0]) & (k2 <= BAND_HI[1])
            out["E2lo"][i, j] = E2[mlo].mean()
            out["Eplo"][i, j] = Ep[mlo].mean()
            out["E2hi"][i, j] = E2[mhi].mean()
            out["Ephi"][i, j] = Ep[mhi].mean()
            R = [trapezoid(t[c], k2) for c in ("E1", "E2", "E3")]
            out["u2frac"][i, j] = R[1] / sum(R)
    fn = os.path.join(os.path.dirname(__file__),
                      "bands_" + os.path.basename(root.rstrip("/")) + ".npz")
    np.savez_compressed(fn, strains=np.array(STRAINS),
                        names=np.array([os.path.basename(d) for d in dirs]),
                        **out)
    print("wrote", fn, f"({n} rlz)")
