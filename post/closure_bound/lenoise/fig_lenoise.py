#!/usr/bin/env python3
"""Part III (main_v3.tex sec:lenoise): the consequence of the computed
distortion for the Amiet leading-edge-noise input.

Chain (Gate-B reduction, geometry-independent): in Amiet's theory the
far-field PSD factorises as [geometry, |L|^2](K_x) * Phi_ww(K_x,0) *
l_y(K_x) with l_y = pi Phi_ww(K_x,0)/Int Phi_ww(K_x,K_y) dK_y.  For a
CHANGE of inflow spectrum at fixed aerofoil/flow/observer, everything
but the turbulence kernel

    K(K_x) = pi * Phi_ww(K_x,0)^2 / Int Phi_ww(K_x,K_y) dK_y

cancels exactly, so
    Delta-SPL(K_x) = 10 log10[ K^distorted / K^frozen ].

Inputs: the RDT-distorted vK representation of the MODEL's strained
line spectra — two fitted parameters (A0, ke) per strain from the
clean 1024-realisation ensembles (rdt_family_fit_table_S1/S20, the
fits behind the distance envelope of sec:envelope; rms residual 13-24%
= the quoted spectral-shape uncertainty).  Phi_ww(K_x,K_y; e) =
A0 * Int Phi_22^{RDT}(K_x/ke, k_2/ke, K_y/ke; e) dk_2 with Phi_22 from
rdt_kernel (exact Cauchy plane-strain distortion of an isotropic vK).
The FROZEN baseline is the same inflow's e=0 fit carried undistorted —
the standard prescribed-spectrum practice.

Axes: K_x in units of the undistorted energy-containing wavenumber
ke(0).  Outputs fig_lenoise.{pdf,png} (canon) + printed band levels.
"""
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import rdt_kernel as rk  # noqa: E402
from figstyle_jfm import COL, FULL, panel, plt, save  # noqa: E402

NK2 = 240      # k2-integral nodes (symmetric log grid)
NKY = 160      # Ky-integral nodes
HALF = 9.0     # log half-width (decades e-folds) for integrals


def read_fits(tag):
    fits = {}
    for line in open(os.path.join(HERE,
                                  f"rdt_family_fit_table_{tag}.txt")):
        m = re.match(r"\s*([\d.]+)\s+([\d.]+)\s+([\d.eE+-]+)\s+[\d.]+",
                     line)
        if m:
            fits[float(m.group(1))] = (float(m.group(2)),
                                       float(m.group(3)))
    return fits


def sym_log_grid(scale, half=HALF, n=NK2):
    """Symmetric grid: +/- log-spaced points around 0 (excludes 0)."""
    x = scale * np.exp(np.linspace(-half, np.log(30.0), n // 2))
    return np.concatenate([-x[::-1], x])


def phi_ww(kx, ky, e, ke, A0):
    """Phi_ww(kx, ky; e) = A0 * Int Phi22^RDT(kx/ke, k2/ke, ky/ke) dk2."""
    k2 = sym_log_grid(ke)
    P22, _ = rk.rdt_phi_components(kx / ke, k2 / ke, ky / ke, e)
    return A0 * np.trapezoid(P22, k2)


def noise_kernel(kx, e, ke, A0):
    p0 = phi_ww(kx, 0.0, e, ke, A0)
    ky = sym_log_grid(ke, n=NKY)
    integ = np.array([phi_ww(kx, kyv, e, ke, A0) for kyv in ky])
    I = np.trapezoid(integ, ky)
    return np.pi * p0 ** 2 / I


def variance(e, ke, A0):
    """Int Int Phi_ww dKx dKy - the level carried by the spectrum; the
    noise kernel scales linearly with it, so dividing it out isolates
    the SHAPE (reallocation) part of Delta-SPL."""
    kx = sym_log_grid(ke, n=NKY)
    ky = sym_log_grid(ke, n=NKY)
    rows = np.array([np.trapezoid(
        np.array([phi_ww(kxv, kyv, e, ke, A0) for kyv in ky]), ky)
        for kxv in kx])
    return np.trapezoid(rows, kx)


def main():
    fig, axs = plt.subplots(1, 2, figsize=(FULL, 2.3), sharey=True)
    for ax, tag, ttl in ((axs[0], "S1",
                          r"$Sk_t/\varepsilon\approx0.4$"),
                         (axs[1], "S20",
                          r"$Sk_t/\varepsilon\approx8$")):
        fits = read_fits(tag)
        ke0, A00 = fits[0.0]
        kx = ke0 * np.geomspace(0.04, 12.0, 34)
        K0 = np.array([noise_kernel(k, 0.0, ke0, A00) for k in kx])
        cols = [COL["odt"], COL["phi33"], COL["phi11"], COL["rdt"]]
        # level term: the MODEL's measured upwash-variance ratio (median
        # line spectra of the clean ensembles), not the fit amplitude
        sp = np.load(os.path.join(HERE, f"spectra_gateA_{tag}.npz"),
                     allow_pickle=True)
        # dumps are ordered e = 0, 0.5, 1, 1.5, 2 in both ensembles
        # (the 'strains' array stores dump TIMES, t = 0.4 + e/S)
        def r22(e):
            i = int(round(e / 0.5))
            k2 = sp[f"k2_e{i}"]
            return np.trapezoid(sp[f"med_E2_e{i}"], k2)
        R220 = r22(0.0)
        V0 = variance(0.0, ke0, A00)
        for (e, (ke, A0)), col in zip(
                [(e, fits[e]) for e in (0.5, 1.0, 1.5, 2.0) if e in fits],
                cols):
            Ke = np.array([noise_kernel(k, e, ke, A0) for k in kx])
            dspl = 10.0 * np.log10(Ke / K0)
            lev = 10.0 * np.log10(r22(e) / R220)
            shape = dspl - lev
            ax.plot(kx / ke0, shape, "-", color=col, ms=2.5, marker="o",
                    label=f"$e={e:.1f}$" if ax is axs[0] else None)
            ax.plot(kx / ke0, dspl, ls=(0, (1, 1.5)), color=col,
                    lw=0.8, alpha=0.6)
            print(f"{tag} e={e}: level {lev:+.1f} dB; shape-only at "
                  "kx/ke0=0.05/0.3/1/3/10: "
                  + " ".join(f"{np.interp(v, kx / ke0, shape):+.1f}"
                             for v in (0.05, 0.3, 1, 3, 10)))
        ax.axhline(0, color="0.6", lw=0.6, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel(r"$K_x/\kappa_e(0)$")
        ax.text(0.5, 0.96, ttl, transform=ax.transAxes, ha="center",
                va="top", fontsize=8)
    axs[0].set_ylabel(r"$\Delta\mathrm{SPL}$ [dB]")
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, ncol=4, loc="lower center",
               bbox_to_anchor=(0.5, -0.14), columnspacing=1.2)
    for j, ax in enumerate(axs):
        panel(ax, "ab"[j], y=0.14)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_lenoise"))


if __name__ == "__main__":
    main()
