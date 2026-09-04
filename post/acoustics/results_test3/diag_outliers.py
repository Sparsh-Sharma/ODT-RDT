#!/usr/bin/env python3
"""Outlier forensics for the Test-3 CARO campaign.

The 1024-rlz table flipped the 64-rlz story (baseline A_high collapsed below 1
with CIs that WIDENED vs 64 rlz).  A = <E2>/<Eperp> is a ratio of ensemble
means, so a few realizations with runaway band energy dominate.  This script
quantifies that: per-realization band means, top offenders by |log energy|
deviation, first-64-only vs full-ensemble A, trimmed A, and median-of-ratios.

    python3 diag_outliers.py <root> [<root2> ...]      # e.g. data/homogeneousStrain2
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import odt_io                                          # noqa: E402

STRAINS = [2.0, 3.9]
BAND_LO = (30.0, 100.0)
BAND_HI = (300.0, 800.0)
NU = 2048


def per_rlz(root, e):
    dirs = sorted(glob.glob(os.path.join(root, "data", "data_*")))
    rows, names = [], []
    for d in dirs:
        t = odt_io.load_fit_target(d, strain=e, Nu=NU)
        k2 = t["k2"]
        Ep = 0.5 * (t["E1"] + t["E3"])
        E2 = t["E2"]
        row = []
        for b in (BAND_LO, BAND_HI):
            m = (k2 >= b[0]) & (k2 <= b[1])
            row += [np.mean(E2[m]), np.mean(Ep[m])]
        rows.append(row)
        names.append(os.path.basename(d))
    return np.array(rows), names


def report(root):
    print(f"=== {root} ===")
    for e in STRAINS:
        r, names = per_rlz(root, e)
        E2lo, Eplo, E2hi, Ephi = r.T
        n = len(names)
        Afull = E2hi.mean() / Ephi.mean()
        A64 = E2hi[:64].mean() / Ephi[:64].mean()
        # per-realization ratio distribution
        Ai = E2hi / Ephi
        med = np.median(Ai)
        # outliers: total high-band energy vs ensemble median
        tot = E2hi + 2.0 * Ephi
        dev = np.log10(tot / np.median(tot))
        order = np.argsort(dev)[::-1]
        # trimmed ratio-of-means (drop top 1% by tot)
        keep = np.sort(order[int(0.01 * n):])
        Atrim = E2hi[keep].mean() / Ephi[keep].mean()
        print(f" e={e:3.1f}  n={n}  A_high: full={Afull:.3f}  first64={A64:.3f}"
              f"  trim1%={Atrim:.3f}  median(E2/Ep)={med:.3f}")
        print(f"   top offenders by high-band energy (log10 dev from median):")
        for i in order[:8]:
            print(f"     {names[i]}  dev={dev[i]:+6.2f}  E2hi={E2hi[i]:.3e}"
                  f"  Ephi={Ephi[i]:.3e}  E2/Ep={Ai[i]:.2f}")
        print(f"   count dev>1 (10x median): {(dev > 1).sum()},"
              f"  dev>2 (100x): {(dev > 2).sum()},  dev>4: {(dev > 4).sum()}")


if __name__ == "__main__":
    roots = sys.argv[1:] or ["../../../data/homogeneousStrain2",
                             "../../../data/homogeneousStrain2A"]
    for root in roots:
        report(root)
