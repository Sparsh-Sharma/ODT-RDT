#!/usr/bin/env python3
"""Manuscript figure fig_cmk_rdt (main_v3.tex sec:cmk): the
rapid-distortion limit against the analytic reference used by
Chen, Meneveau & Katz (2006) — eddy events suppressed, inviscid.
(a) component spectral-centroid migration vs the rigid-translation law,
to the end of CMK's straining phase (D=3, e=2.20);
(b) kinetic-energy amplification vs exact RDT (wavevector-ensemble
integration) and the model closure trajectory.
Data: hsA2_ens.npz (case hsA2fig rerun on caro 2026-09-06)."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "verification"))
from figstyle_jfm import COL, FULL, panel, plt, save  # noqa: E402
from fig_rdt_moments import (exact_rdt, integrate_model,  # noqa: E402
                             strain_tensor)

ECMK = 2.20   # end of CMK's straining phase: D = exp(e/2) = 3


def main():
    d = np.load(os.path.join(HERE, "hsA2_ens.npz"))
    t = d["t"]
    m = t <= ECMK + 1e-9
    A = strain_tensor("plane")
    te, Re = exact_rdt(A, emax=ECMK)
    kte = 0.5 * (Re[0, 0] + Re[1, 1] + Re[2, 2])
    tm, Rm = integrate_model(A, "LRR", emax=ECMK)
    ktm = 0.5 * (Rm[0, 0] + Rm[1, 1] + Rm[2, 2])

    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.3))
    ee = np.linspace(0, ECMK, 100)
    a.plot(ee, np.exp(0.5 * ee), "k--", lw=0.9,
           label=r"rigid translation $\mathrm{e}^{-A_{22}e}$")
    a.plot(t[m], d["c1_mean"][m] / d["c1_mean"][0], "o",
           color=COL["phi11"], ms=2.8, mfc="none",
           label=r"$\bar\kappa_1$")
    a.plot(t[m], d["c2_mean"][m] / d["c2_mean"][0], "s",
           color=COL["phi22"], ms=2.8, mfc="none",
           label=r"$\bar\kappa_2$")
    a.axhline(3.0, color="0.7", lw=0.5, ls=":")
    a.set_ylim(0.95, 3.3)
    a.text(1.9, 3.05, r"$D=3$", fontsize=7, color="0.4")
    a.set_xlabel(r"total strain $e = S\,t$")
    a.set_ylabel(r"$\bar\kappa(e)/\bar\kappa(0)$")
    a.legend(loc="upper left")
    b.plot(te, kte / kte[0], color=COL["rdt"], lw=1.2,
           label="exact RDT")
    b.plot(tm, ktm / ktm[0], color="k", ls=(0, (6, 2)), lw=0.9,
           label="LRR closure")
    b.plot(t[m], d["kt_mean"][m] / d["kt_mean"][0], "o",
           color=COL["odt"], ms=2.8, mfc="none", label="ODT solver")
    b.set_xlabel(r"total strain $e = S\,t$")
    b.set_ylabel(r"$k_t(e)/k_t(0)$")
    b.legend(loc="upper left")
    for j, ax in enumerate((a, b)):
        ax.set_xlim(0, ECMK)
        panel(ax, "ab"[j], x=0.88, y=0.15)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_cmk_rdt"))
    print("centroid at e=2.2:",
          np.interp(ECMK, t, d["c2_mean"] / d["c2_mean"][0]))
    print("kt run/LRR/exact at e=2.2: "
          f"{np.interp(ECMK, t, d['kt_mean'] / d['kt_mean'][0]):.3f} "
          f"{np.interp(ECMK, tm, ktm / ktm[0]):.3f} "
          f"{np.interp(ECMK, te, kte / kte[0]):.3f}")


if __name__ == "__main__":
    main()
