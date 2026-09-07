#!/usr/bin/env python3
"""Manuscript figure fig_clock (Part II, sec:clock): the independent
relaxation clock against the eddy-locked family on the transmission
diagnostic.  (a) fine-scale-band imbalance A_hi(e): the eddy-locked
variants fade beyond e~3, the clock holds at every strain; the uncapped
clock included to show the overcorrection.  (b) the price panel:
upwash energy fraction u2^2/2kt(e) — the uncapped clock crushes the
physical one-point anisotropy, the size cap restores it.
Data: newbands/*.npz (caro campaigns 2026-09-07/08, 1024 paired-seed
realisations each, one binary)."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from figstyle_jfm import BLACK, GOLD, GREEN, PURPLE, RED, FULL, panel, plt, save  # noqa: E402,E501

NB = os.path.join(HERE, "newbands")
CASES = [
    ("standard ODT", "bands_homogeneousStrain2.npz", BLACK, "o", "-"),
    ("depth tied to eddy size", "bands_homogeneousStrain2KS.npz",
     PURPLE, "v", "-"),
    ("clock, uncapped", "bands_homogeneousStrain2RC1.npz", GOLD, "s",
     (0, (4, 2))),
    ("clock, size-capped", "bands_homogeneousStrain2RCS1.npz", RED,
     "^", "-"),
    ("clock, size-capped, $4\\times$", "bands_homogeneousStrain2RCS4.npz",
     GREEN, "d", (0, (1, 1.5))),
]
NBOOT = 2000


def med_ci(x, rng):
    m = np.median(x)
    idx = rng.integers(0, x.size, (NBOOT, x.size))
    mb = np.median(x[idx], axis=1)
    return m, np.percentile(mb, 2.5), np.percentile(mb, 97.5)


def main():
    rng = np.random.default_rng(0)
    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.3), sharex=True)
    for lab, fn, col, mk, ls in CASES:
        d = np.load(os.path.join(NB, fn))
        es = d["strains"]
        Ahi = d["E2hi"] / d["Ephi"]
        m, lo, hi = np.array([med_ci(Ahi[:, j], rng)
                              for j in range(len(es))]).T
        a.errorbar(es, m, yerr=[m - lo, hi - m], fmt=mk, ls=ls,
                   color=col, capsize=1.8, capthick=0.6, elinewidth=0.6,
                   ms=3, lw=0.9, label=lab)
        u2 = np.median(d["u2frac"], axis=0)
        b.plot(es, u2, marker=mk, ls=ls, color=col, ms=3, lw=0.9)
    a.set_ylabel(r"$\mathcal{A}_{\rm hi}$")
    a.set_xlabel(r"total strain $e$")
    b.set_ylabel(r"$u_2^2/2k_t$")
    b.set_xlabel(r"total strain $e$")
    h, lab = a.get_legend_handles_labels()
    fig.legend(h, lab, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, -0.16), columnspacing=1.0)
    panel(a, "a")
    panel(b, "b")
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_clock"))


if __name__ == "__main__":
    main()
