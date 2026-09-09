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
from figstyle_jfm import FULL, plt, save  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

x = np.linspace(0.0, 1.0, 400)
# monochrome: kernel -> marker; coefficient -> line style
# (a solid black, c dashed grey).
MK = {r"$g_1$": "o", r"$g_2$": "s", r"$g_3$": "^"}
GA = {r"$g_1$": 0.125 + 0.75 * x - 0.875 * x**2,
      r"$g_2$": -1.5 * x * (1 - x),
      r"$g_3$": -0.125 + 0.75 * x - 0.625 * x**2}
GC = {r"$g_1$": 0.875 * x * (1 - x)**2,
      r"$g_2$": -1.5 * x * (1 - x)**2,
      r"$g_3$": 0.625 * x * (1 - x)**2}


def main():
    fig, ax = plt.subplots(figsize=(0.72 * FULL, 2.5))
    for i, lab in enumerate((r"$g_1$", r"$g_2$", r"$g_3$")):
        ev = (40 * i + 20, 120)
        ax.plot(x, GA[lab], color="k", ls="-", lw=1.2, marker=MK[lab],
                ms=3.4, mfc="k", mew=0.6, markevery=ev)
        ax.plot(x, GC[lab], color="0.5", ls=(0, (5, 2)), lw=0.9,
                marker=MK[lab], ms=3.0, mfc="none", mew=0.7, markevery=ev)
    ax.axhline(0, color="0.6", lw=0.6, ls=":")
    ax.set_xlabel(r"$x=\kappa_2^2/\kappa^2$")
    ax.set_ylabel(r"kernel value")
    kern_h = [Line2D([], [], color="k", ls="-", lw=1.2, marker=MK[k],
                     ms=3.4, mfc="k", mew=0.6, label=k)
              for k in (r"$g_1$", r"$g_2$", r"$g_3$")]
    coef_h = [Line2D([], [], color="k", ls="-", lw=1.2,
                     label=r"$a$ coefficient"),
              Line2D([], [], color="0.5", ls=(0, (5, 2)), lw=0.9,
                     label=r"$c$ coefficient")]
    ax.legend(handles=kern_h + coef_h, ncol=2, loc="lower left",
              columnspacing=1.0, fontsize=6.8)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_pi_kernels"))


if __name__ == "__main__":
    main()
