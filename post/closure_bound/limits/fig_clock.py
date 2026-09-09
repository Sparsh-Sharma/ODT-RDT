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
from figstyle_jfm import FULL, panel, plt, save  # noqa: E402

NB = os.path.join(HERE, "newbands")
# monochrome: each intervention -> distinct line style + marker (all black)
CASES = [
    ("SC-ODT (baseline)", "bands_homogeneousStrain2.npz", "o", "-"),
    ("depth tied to eddy size", "bands_homogeneousStrain2KS.npz",
     "v", (0, (6, 2))),
    ("clock, uncapped", "bands_homogeneousStrain2RC1.npz", "s",
     (0, (1, 1.4))),
    ("clock-relaxed ODT", "bands_homogeneousStrain2RCS1.npz",
     "^", (0, (5, 1.6, 1, 1.6))),
    ("clock-relaxed ODT, $4\\times$", "bands_homogeneousStrain2RCS4.npz",
     "D", (0, (4, 1.4, 1, 1.4, 1, 1.4))),
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
    for lab, fn, mk, ls in CASES:
        d = np.load(os.path.join(NB, fn))
        es = d["strains"]
        Ahi = d["E2hi"] / d["Ephi"]
        m, lo, hi = np.array([med_ci(Ahi[:, j], rng)
                              for j in range(len(es))]).T
        a.errorbar(es, m, yerr=[m - lo, hi - m], fmt=mk, ls=ls,
                   color="k", capsize=1.6, capthick=0.5, elinewidth=0.5,
                   ms=3, mfc="k", lw=0.9, label=lab)
        u2 = np.median(d["u2frac"], axis=0)
        b.plot(es, u2, marker=mk, ls=ls, color="k", ms=3, mfc="k", lw=0.9)
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
