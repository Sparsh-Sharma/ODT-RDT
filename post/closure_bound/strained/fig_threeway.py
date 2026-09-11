#!/usr/bin/env python3
"""Manuscript figure fig_threeway (main_v3.tex sec:threeway), plot-only:
reads threeway.npz written by threeway.py (which does the computation)
and renders the three-way comparison in the JFM house style at true
\\textwidth.  Panels: (a) strain-induced upwash anisotropy Delta b_22;
(b) Delta b_11 (solid) and Delta b_33 (dotted); (c) transverse splitting.
e=1 solid, e=0.5 dashed (RDT and DNS only)."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from figstyle_jfm import FULL, panel, plt, save  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

D = np.load(os.path.join(HERE, "threeway.npz"), allow_pickle=True)
EDGES = D["edges"]
XC = np.sqrt(EDGES[:-1] * EDGES[1:])

# monochrome: system -> marker, strain -> solid(e=1)/dashed(e=0.5),
# db33 -> dotted.  All black; SC-ODT is the model (was "ISO").
MK = {"RDT": "D", "DNS": "o", "ISO": "s"}
LBL = {"RDT": "exact RDT, projected", "DNS": r"DNS $128^3$",
       "ISO": "SC-ODT"}
LW = {"RDT": 1.0, "DNS": 1.2, "ISO": 1.0}


def q(s, e, name):
    return D[f"{s}_e{e}_{name}"]


def main():
    # one component per panel (2x2); three systems at the operating strain
    # e=1 only -- one filled marker per system, single line style, no
    # second encoding.  (Strain-rapidity dependence is the subject of the
    # following envelope figure, sec:envelope.)
    # three diagonal components in a row (SC-ODT is not the final model,
    # so its evaluation need not be exhaustive -- Alan, 2026-09-11): the
    # transverse splitting phi11/phi33-1 is derivable from b11, b33 and is
    # described in the text rather than shown.
    fig, axs = plt.subplots(1, 3, figsize=(FULL, 2.3), sharex=True)
    panels = ((axs[0], "db22", r"$\Delta b_{22}(\kappa_2)$"),
              (axs[1], "db11", r"$\Delta b_{11}(\kappa_2)$"),
              (axs[2], "db33", r"$\Delta b_{33}(\kappa_2)$"))
    for ax, name, yl in panels:
        for s in ("RDT", "DNS", "ISO"):
            ax.plot(XC, q(s, 1.0, name), color="k", ls="-", lw=LW[s],
                    marker=MK[s], ms=3.2, mfc="k", mew=0.6)
        ax.axhline(0, color="0.75", lw=0.6, ls=":")
        ax.set_xscale("log")
        ax.set_ylabel(yl)
        ax.set_xlabel(r"$\kappa_2(e)/\kappa_c(0)$")
    for j, ax in enumerate(axs):
        panel(ax, "abc"[j], y=0.15)
    sys_h = [Line2D([], [], color="k", ls="-", lw=1.1, marker=MK[s],
                    ms=3.2, mfc="k", mew=0.6, label=LBL[s])
             for s in ("RDT", "DNS", "ISO")]
    fig.legend(handles=sys_h, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, -0.04), columnspacing=1.6,
               handlelength=2.2, fontsize=7.5)
    fig.tight_layout(pad=0.5, rect=(0, 0.08, 1, 1))
    save(fig, os.path.join(HERE, "fig_threeway"))


if __name__ == "__main__":
    main()
