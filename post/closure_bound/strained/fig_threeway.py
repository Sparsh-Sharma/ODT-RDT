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
from figstyle_jfm import FULL, plt, save, strip  # noqa: E402

D = np.load(os.path.join(HERE, "threeway.npz"), allow_pickle=True)
EDGES = D["edges"]
XC = np.sqrt(EDGES[:-1] * EDGES[1:])

COL = {"RDT": "#8e44ad", "DNS": "k", "ISO": "#2a78d6", "TYPESw": "#eda100"}
LBL = {"RDT": "exact linear RDT, projected",
       "DNS": r"DNS $128^3$",
       "ISO": r"ODT, isotropic kernel",
       "TYPESw": r"ODT, eddy types $p_3{=}0.5$"}
LW = {"RDT": 1.3, "DNS": 1.5, "ISO": 1.2, "TYPESw": 0.8}
AL = {"RDT": 1, "DNS": 1, "ISO": 1, "TYPESw": 0.5}
MK = {"RDT": "D", "DNS": "o", "ISO": "s", "TYPESw": "^"}


def q(s, e, name):
    return D[f"{s}_e{e}_{name}"]


def main():
    fig, axs = plt.subplots(1, 3, figsize=(FULL, 2.1))
    panels = ((axs[0], "db22", r"$\Delta b_{22}(\kappa_2)$"),
              (axs[1], "db11", r"$\Delta b_{11},\ \Delta b_{33}$"),
              (axs[2], "split", r"$\phi_{11}/\phi_{33}-1$"))
    for ax, name, yl in panels:
        for s in ("RDT", "DNS", "ISO", "TYPESw"):
            ax.plot(XC, q(s, 1.0, name), "-", color=COL[s], lw=LW[s],
                    alpha=AL[s], marker=MK[s], ms=2.6,
                    label=(LBL[s] + ", $e{=}1$") if ax is axs[0] else None)
            if s in ("RDT", "DNS"):
                ax.plot(XC, q(s, 0.5, name), "--", color=COL[s],
                        lw=LW[s] * 0.55, alpha=0.8,
                        label=(LBL[s].split(",")[0] + ", $e{=}0.5$")
                        if ax is axs[0] else None)
            if name == "db11":
                ax.plot(XC, q(s, 1.0, "db33"), ":", color=COL[s],
                        lw=LW[s], alpha=AL[s], marker=MK[s], ms=2.2)
        ax.axhline(0, color="0.6", lw=0.6, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel(r"$\kappa_2(e)/\kappa_c(0)$")
        ax.set_ylabel(yl)
        strip(ax)
    for j, ax in enumerate(axs):
        ax.text(0.03, 0.97, "(" + "abc"[j] + ")",
                transform=ax.transAxes, va="top")
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, -0.16), handlelength=1.6,
               columnspacing=1.0)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_threeway"))


if __name__ == "__main__":
    main()
