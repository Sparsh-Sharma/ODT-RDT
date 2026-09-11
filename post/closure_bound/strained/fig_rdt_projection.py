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
from figstyle_jfm import FULL, panel, plt, save  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

# DNS (strained box, 128^3): binned splitting from rdt_projection.npz.
# Drop the first (largest-scale) bin: it holds a single integer kappa_2
# line and its e=1/e=0.5 entries are one-sided.
D = np.load(os.path.join(HERE, "rdt_projection.npz"), allow_pickle=True)
EDGES = D["edges"]
SL = slice(1, None)
XC = np.sqrt(EDGES[:-1] * EDGES[1:])[SL]


def g(key):
    return D[key][SL]


# exact RDT: closed-form Cauchy solution integrated analytically over the
# (k1, k3) plane (rdt_projection_analytic.py), vK initial spectrum
# band-limited to the DNS box -- replaces the 4-seed companions, which
# were a noisy estimate of this integral at the largest scales.
A = np.load(os.path.join(HERE, "rdt_projection_analytic.npz"),
            allow_pickle=True)
# plot over the DNS-resolved range only: beyond x ~ 10 the curves feel the
# box's dealiasing edge, which is a finite-box effect, not linear theory
MA = A["x"] <= 10.0
XA = A["x"][MA]


def h(key):
    return A[key][MA]

# monochrome: component -> marker, strain -> line style + marker fill
MK = {"phi22": "o", "phi11": "s", "phi33": "^"}
LAB = {"phi22": r"$\phi_{22}$", "phi11": r"$\phi_{11}$",
       "phi33": r"$\phi_{33}$"}


def main():
    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.2))

    # ---- (a) shape ratio: component -> marker, strain -> solid/dashed ---
    for ci, name in enumerate(("phi22", "phi11", "phi33")):
        mk = MK[name]
        a.plot(XA, h(f"R_{name}_e1"), color="k", ls="-", lw=1.0,
               marker=mk, ms=3.4, mfc="k", mew=0.8,
               markevery=(3 * ci + 1, 9))
        a.plot(XA, h(f"R_{name}_e0.5"), color="k", ls=(0, (5, 2)),
               lw=0.8, marker=mk, ms=3.4, mfc="none", mew=0.8,
               markevery=(3 * ci + 1, 9))
    a.axhline(1, color="0.5", lw=0.7, ls=":")
    a.set_xscale("log")
    a.set_ylim(0.25, 1.32)
    a.set_xlabel(r"$\kappa_2(e)/\kappa_c(0)$")
    a.set_ylabel(r"shape ratio $R_{nn}$")
    comp_h = [Line2D([], [], color="k", ls="none", marker=MK[n], ms=3.4,
                     mfc="k", mew=0.8, label=LAB[n])
              for n in ("phi22", "phi11", "phi33")]
    strain_h = [
        Line2D([], [], color="k", ls="-", lw=1.0, marker="o", ms=3.4,
               mfc="k", mew=0.8, label=r"$e{=}1$"),
        Line2D([], [], color="k", ls=(0, (5, 2)), lw=0.8, marker="o",
               ms=3.4, mfc="none", mew=0.8, label=r"$e{=}0.5$"),
        Line2D([], [], color="0.5", ls=":", lw=0.7,
               label="rigid translation")]
    a.legend(handles=comp_h + strain_h, ncol=2, columnspacing=0.8,
             loc="lower left", bbox_to_anchor=(0.01, 0.01), fontsize=6.5)

    # ---- (b) transverse splitting: source -> line style, strain -> lw ---
    b.plot(XC, g("split_dns_e1"), color="k", ls="-", lw=1.3,
           marker="s", ms=3.2, mfc="none", mew=0.8, markevery=1)
    b.plot(XA, h("split_e1"), color="k", ls=(0, (5, 2)), lw=1.3)
    b.plot(XC, g("split_dns_e0.5"), color="k", ls="-", lw=0.7,
           marker="s", ms=2.6, mfc="none", mew=0.6, markevery=1)
    b.plot(XA, h("split_e0.5"), color="k", ls=(0, (5, 2)), lw=0.7)
    b.axhline(0, color="0.6", lw=0.6, ls=":")
    b.set_xscale("log")
    b.set_xlabel(r"$\kappa_2(e)/\kappa_c(0)$")
    b.set_ylabel(r"$\phi_{11}/\phi_{33}-1$")
    leg = b.legend(handles=[
        Line2D([], [], color="k", ls="-", lw=1.3, marker="s", ms=3.2,
               mfc="none", mew=0.8, label=r"DNS, $e{=}1$"),
        Line2D([], [], color="k", ls=(0, (5, 2)), lw=1.3,
               label=r"exact RDT, $e{=}1$"),
        Line2D([], [], color="k", ls="-", lw=0.7, marker="s", ms=2.6,
               mfc="none", mew=0.6, label=r"DNS, $e{=}0.5$"),
        Line2D([], [], color="k", ls=(0, (5, 2)), lw=0.7,
               label=r"exact RDT, $e{=}0.5$")],
        loc="upper left", fontsize=6.5, ncol=1,
        frameon=True, framealpha=0.92, edgecolor="none")
    leg.get_frame().set_facecolor("white")
    panel(a, "a", x=0.02, y=0.97)
    panel(b, "b", x=0.02, y=0.97)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_rdt_projection"))


if __name__ == "__main__":
    main()
