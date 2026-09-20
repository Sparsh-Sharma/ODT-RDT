#!/usr/bin/env python3
"""The SC-ODT vs clock-relaxed difference, made visible for a general
reader. A single realisation cannot separate the two (they diverge into
different realisations); the difference is statistical. Over the
1024-realisation ensembles it is clear: the fine-scale-band imbalance
A_hi = E2_hi / Ep_hi (upwash over spanwise energy in the small-scale
band) grows and persists under strain for SC-ODT, whereas the
clock-relaxed model relaxes it back towards balance (A_hi -> 1) at every
strain.

  python fig_hero_diff.py            colour (proposals)
  python fig_hero_diff.py --mono     greyscale (paper)
Data: post/closure_bound/limits/newbands/bands_homogeneousStrain2*.npz
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from figstyle_jfm import FULL, plt, save, COL  # noqa: E402

MONO = "--mono" in sys.argv
NB = os.path.join(HERE, "limits", "newbands")
NBOOT = 2000
CASES = [("SC-ODT", "bands_homogeneousStrain2.npz",
          "0.35" if MONO else COL["odt"], "o", (0, (5, 2))),
         ("clock-relaxed ODT", "bands_homogeneousStrain2RCS1.npz",
          "k" if MONO else COL["b22"], "s", "-")]


def med_ci(x, rng):
    m = np.median(x)
    idx = rng.integers(0, x.size, (NBOOT, x.size))
    mb = np.median(x[idx], axis=1)
    return m, np.percentile(mb, 2.5), np.percentile(mb, 97.5)


def main():
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(0.66 * FULL, 2.7))
    ax.axhline(1.0, color="0.7", lw=0.6, ls=":")
    ax.text(0.06, 1.02, "balanced (isotropic) fine scales",
            fontsize=6.4, color="0.45", va="bottom")
    for lab, fn, col, mk, ls in CASES:
        d = np.load(os.path.join(NB, fn))
        es = d["strains"]
        Ahi = d["E2hi"] / d["Ephi"]
        m, lo, hi = np.array([med_ci(Ahi[:, j], rng)
                              for j in range(len(es))]).T
        ax.fill_between(es, lo, hi, color=col, alpha=0.15, lw=0)
        ax.plot(es, m, color=col, ls=ls, marker=mk, ms=4, mfc=col,
                lw=1.6, label=lab)
    ax.annotate("SC-ODT: the strain-driven\nimbalance persists",
                xy=(3.2, 1.82), xytext=(1.75, 1.62), fontsize=7,
                color=(("0.3") if MONO else COL["odt"]), ha="left",
                arrowprops=dict(arrowstyle="-|>",
                                color=("0.4" if MONO else COL["odt"]),
                                lw=0.7))
    ax.annotate("clock-relaxed: the clock\nrelaxes it back at every strain",
                xy=(3.0, 1.36), xytext=(1.15, 1.14), fontsize=7,
                color=("k" if MONO else COL["b22"]), ha="left",
                arrowprops=dict(arrowstyle="-|>",
                                color=("0.2" if MONO else COL["b22"]),
                                lw=0.7))
    ax.set_xlabel(r"total strain $e = St$")
    ax.set_ylabel(r"fine-scale imbalance $\mathcal{A}_{\rm hi}$")
    ax.set_xlim(-0.1, 4.1)
    ax.set_ylim(0.95, 2.15)
    ax.legend(loc="upper left", fontsize=7)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_hero_diff" + ("_mono" if MONO else "")))


if __name__ == "__main__":
    main()
