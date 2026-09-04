#!/usr/bin/env python3
"""Robust (outlier-immune) baseline vs Option-A comparison from bands_*.npz.

The per-realization band energies are heavy-tailed (single realizations carry
100-2000x the median band energy), so the ratio-of-ensemble-means A used in
compare_optionA.py does not converge even at 1024 realizations.  Here:

  A_med  = median_i( E2_i / Eperp_i )          per band, bootstrap CI of median
  T_med  = median_i( A_high,i / A_low,i )      per-realization transmission
  paired = median_i( log A_high,i^optA - log A_high,i^base )  same-seed pairs

    python3 robust_optionA.py bands_homogeneousStrain2.npz bands_homogeneousStrain2A.npz
"""
import os
import sys

import numpy as np

NBOOT = 4000


def med_ci(x, rng):
    m = np.median(x)
    idx = rng.integers(0, x.size, (NBOOT, x.size))
    mb = np.median(x[idx], axis=1)
    return m, np.percentile(mb, 2.5), np.percentile(mb, 97.5)


def fmt(ci):
    return f"{ci[0]:5.2f} [{ci[1]:4.2f},{ci[2]:4.2f}]"


def main(f_base, f_opta):
    here = os.path.dirname(os.path.abspath(__file__))
    b = np.load(f_base, allow_pickle=True)
    a = np.load(f_opta, allow_pickle=True)
    strains = b["strains"]
    rng = np.random.default_rng(0)

    lines = [
        f"baseline: {os.path.basename(f_base)} ({b['E2hi'].shape[0]} rlz)   "
        f"option A: {os.path.basename(f_opta)} ({a['E2hi'].shape[0]} rlz)",
        "robust statistics: median over realizations, 95% bootstrap CI of the median",
        "",
        f"{'':>5} | {'baseline':^42} | {'option A':^42}",
        f"{'e':>5} | {'A_low':^17} {'A_high':^17} {'u2/2kt':>6} |"
        f" {'A_low':^17} {'A_high':^17} {'u2/2kt':>6}",
    ]
    res = {}
    for c, d in (("base", b), ("optA", a)):
        res[c] = {
            "Alo": d["E2lo"] / d["Eplo"],       # (n, ne)
            "Ahi": d["E2hi"] / d["Ephi"],
            "u2": d["u2frac"],
        }
    for j, e in enumerate(strains):
        row = f"{e:5.1f} |"
        for c in ("base", "optA"):
            r = res[c]
            row += (f" {fmt(med_ci(r['Alo'][:, j], rng))}"
                    f" {fmt(med_ci(r['Ahi'][:, j], rng))}"
                    f" {np.median(r['u2'][:, j]):6.3f}")
            if c == "base":
                row += " |"
        lines.append(row)

    lines += ["", "downscale transmission, median_i(A_high,i/A_low,i) [95% CI]:"]
    for j, e in enumerate(strains):
        if e == 0:
            continue
        tb = med_ci(res["base"]["Ahi"][:, j] / res["base"]["Alo"][:, j], rng)
        ta = med_ci(res["optA"]["Ahi"][:, j] / res["optA"]["Alo"][:, j], rng)
        lines.append(f"  e={e:3.1f}:  baseline {fmt(tb)}   option A {fmt(ta)}")

    lines += ["", "paired same-seed contrast, median_i[log10(A^optA/A^base)] [95% CI]:",
              "  (negative = option A less anisotropic than baseline for the same IC)"]
    nb = min(res["base"]["Ahi"].shape[0], res["optA"]["Ahi"].shape[0])
    for j, e in enumerate(strains):
        if e == 0:
            continue
        out = f"  e={e:3.1f}: "
        for band in ("Alo", "Ahi"):
            d = np.log10(res["optA"][band][:nb, j] / res["base"][band][:nb, j])
            m, lo, hi = med_ci(d, rng)
            star = " *" if (lo > 0 or hi < 0) else "  "
            out += f"  {band[1:]}: {m:+6.3f} [{lo:+6.3f},{hi:+6.3f}]{star}"
        lines.append(out)

    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(here, "optionA_robust_table.txt"), "w").write(txt + "\n")


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    fb = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "bands_homogeneousStrain2.npz")
    fa = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "bands_homogeneousStrain2A.npz")
    main(fb, fa)
