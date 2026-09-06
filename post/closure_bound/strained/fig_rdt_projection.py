#!/usr/bin/env python3
"""Manuscript figure fig_rdt_projection (main_v3.tex sec:spec-linref),
plot-only: reads rdt_projection.npz written by rdt_projection.py and
renders the two panels in the JFM house style at true \\textwidth.
(a) shape ratio R_nn of the exactly distorted, line-projected component
spectra to a rigid translation; (b) transverse splitting vs wavenumber,
exact RDT against the strained-box DNS."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from figstyle_jfm import COL, FULL, panel, plt, save  # noqa: E402

D = np.load(os.path.join(HERE, "rdt_projection.npz"), allow_pickle=True)
EDGES = D["edges"]
XC = np.sqrt(EDGES[:-1] * EDGES[1:])


def main():
    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.2))
    for name in ("phi22", "phi11", "phi33"):
        col = COL[name]
        a.plot(XC, D[f"R_{name}_e1"], "o-", color=col,
               label=r"$\phi_{%s}$, $e{=}1$" % name[3:])
        a.plot(XC, D[f"R_{name}_e0.5"], "o--", color=col, ms=2.4,
               lw=0.8, alpha=0.6,
               label=r"$\phi_{%s}$, $e{=}0.5$" % name[3:])
    a.axhline(1, color="k", lw=0.7, ls=":", label="rigid translation")
    a.set_xscale("log")
    a.set_xlabel(r"$\kappa_2(e)/\kappa_c(0)$")
    a.set_ylabel(r"shape ratio $R_{nn}$")
    a.legend(ncol=2, columnspacing=0.8, loc="upper right")
    b.plot(XC, D["split_dns_e1"], "-", color=COL["dns"], lw=1.4,
           label=r"DNS, $e{=}1$")
    b.plot(XC, D["split_rdt_e1"], "-", color=COL["rdt"], lw=1.3,
           label="exact RDT, $e{=}1$")
    b.plot(XC, D["split_dns_e0.5"], "--", color=COL["dns"], lw=0.8,
           label=r"DNS, $e{=}0.5$")
    b.plot(XC, D["split_rdt_e0.5"], "--", color=COL["rdt"], lw=0.8,
           label="exact RDT, $e{=}0.5$")
    b.axhline(0, color="0.75", lw=0.6, ls=":")
    b.set_xscale("log")
    b.set_xlabel(r"$\kappa_2(e)/\kappa_c(0)$")
    b.set_ylabel(r"$\phi_{11}/\phi_{33}-1$")
    b.legend(loc="lower left")
    for j, ax in enumerate((a, b)):
        panel(ax, "ab"[j], y=0.08)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_rdt_projection"))


if __name__ == "__main__":
    main()
