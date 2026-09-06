#!/usr/bin/env python3
"""Manuscript figure fig_pi_kernels (main_v3.tex sec:lrr-reduction):
the rapid pressure--strain closure kernels g_1, g_2, g_3 of
eqs. (g1)-(g3) against x = kappa_2^2/kappa^2.  Solid: coefficient of the
isotropic-like scalar a; dashed faint: coefficient of the anisotropy
scalar c.  Verified symbolically against A_kl M_nnkl (2026-09-06).
Output at 0.72\\textwidth."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from figstyle_jfm import FULL, plt, save, strip  # noqa: E402

x = np.linspace(0.0, 1.0, 400)
GA = {r"$g_1$": (0.125 + 0.75 * x - 0.875 * x**2, "#eb6834"),
      r"$g_2$": (-1.5 * x * (1 - x), "#2a78d6"),
      r"$g_3$": (-0.125 + 0.75 * x - 0.625 * x**2, "#1baf7a")}
GC = {r"$g_1$": 0.875 * x * (1 - x)**2,
      r"$g_2$": -1.5 * x * (1 - x)**2,
      r"$g_3$": 0.625 * x * (1 - x)**2}


def main():
    fig, ax = plt.subplots(figsize=(0.72 * FULL, 2.5))
    for lab, (ga, col) in GA.items():
        ax.plot(x, ga, "-", color=col, label=lab + " ($a$ coefficient)")
        ax.plot(x, GC[lab], "--", color=col, alpha=0.45, lw=0.9,
                label=lab + " ($c$ coefficient)")
    ax.axhline(0, color="0.6", lw=0.6, ls=":")
    ax.set_xlabel(r"$x=\kappa_2^2/\kappa^2$")
    ax.set_ylabel(r"kernel value")
    ax.legend(ncol=2, loc="lower left", handlelength=1.8,
              labelspacing=0.25, columnspacing=1.0)
    strip(ax)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_pi_kernels"))


if __name__ == "__main__":
    main()
