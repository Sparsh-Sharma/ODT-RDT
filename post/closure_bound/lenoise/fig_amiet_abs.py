#!/usr/bin/env python3
"""Absolute Amiet leading-edge-noise prediction: standard practice
(frozen von Karman inflow spectrum) against the model's computed
distorted spectrum, at the paper's reference configuration.

Chain: compact lift-dipole far-field PSD with the exact Sears response
    S_pp(omega) = (omega x3/(4 pi c0 sigma^2))^2 (pi rho0 U c)^2 |S|^2
                  * (Phi_ww(omega/U,0)/U) * d * ell_y ,
the normalisation-audited form of the author's len_vsdb model
(amiet_reference.far_field_psd_dipole; the 2*pi^2 audit) with the exact
Sears function (gust_response.sears, Theodorsen/Hankel) — both
reproduced here verbatim to keep this script self-contained.

Inflow spectra: the RDT-distorted vK representation of the model's
strained ensembles (as in fig_lenoise; Sk/eps~0.4 fits), bridged to
physical units by (i) the vK relation kappa_e = 0.7468/Lambda fixing
the length unit and (ii) the upwash variance u'^2 fixing the amplitude
at e=0; the strain-level change uses the measured line-spectrum
variance ratio.  The e=0 member IS the isotropic von Karman spectrum,
so the frozen baseline is exactly "Amiet + vK" standard practice.

Configuration (paper sec:spec-op): U=20 m/s, u'/U=6%, Lambda=55 mm,
chord 0.386 m (Re_c=5.1e5), span 0.8 m, observer 1 m overhead,
M=U/c0.  Output: SPL(f) dB/Hz re 20 uPa, fig_amiet_abs.{pdf,png}.
"""
import os
import sys

import numpy as np
from scipy.special import jv, yv

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import fig_lenoise as fl  # noqa: E402  (phi_ww, read_fits, sym_log_grid)
from figstyle_jfm import COL, FULL, plt, save  # noqa: E402

RHO0, C0, PREF = 1.225, 340.0, 20e-6      # len_vsdb constants
U = 20.0
UP = 0.06 * U                              # u' (upwash rms)
LAMBDA = 0.055                             # integral scale, m
CHORD = 0.386                              # Re_c = 5.1e5 at U=20
SPAN = 0.8
XOBS = (0.0, 0.0, 1.0)                     # 1 m overhead
M = U / C0
CVK = 0.746834                             # kappa_e = CVK / Lambda (vK)


def sears(k):
    """Exact Sears function (len_vsdb gust_response, Theodorsen/Hankel)."""
    k = np.asarray(k, float)
    H1 = jv(1, k) - 1j * yv(1, k)
    H0 = jv(0, k) - 1j * yv(0, k)
    C = H1 / (H1 + 1j * H0)
    return C * (jv(0, k) - 1j * jv(1, k)) + 1j * jv(1, k)


def far_field_psd_dipole(x, omega, S, phi_ww, chord, U, span, M, ell_y):
    """len_vsdb amiet_reference.far_field_psd_dipole, verbatim."""
    x1, x2, x3 = x
    sigma2 = x1 ** 2 + (1 - M ** 2) * (x2 ** 2 + x3 ** 2)
    dipole2 = (omega * x3 / (4 * np.pi * C0 * sigma2)) ** 2
    lift_per_w2 = (np.pi * RHO0 * U * chord) ** 2 * np.abs(S) ** 2
    S_ww = phi_ww / U
    return dipole2 * lift_per_w2 * S_ww * span * ell_y


def vk_phi1d(Kx, w2, Lam):
    """Analytic 1-D von Karman upwash spectrum, two-sided (Amiet 1975)."""
    kh = np.asarray(Kx, float) * Lam / CVK
    return (w2 * Lam / (6 * np.pi)) * (3 + 8 * kh ** 2) \
        / (1 + kh ** 2) ** (11.0 / 6.0)


def liepmann_phi1d(Kx, w2, Lam):
    """len_vsdb turbulence.liepmann_phi_ww, verbatim (two-sided)."""
    x = (Lam * np.asarray(Kx, float)) ** 2
    return (w2 * Lam / (2 * np.pi)) * (1 + 3 * x) / (1 + x) ** 2


def main():
    fits = fl.read_fits("S1")
    ke0, A00 = fits[0.0]
    # unit bridge: model length unit in metres
    ke_phys = CVK / LAMBDA                  # 13.58 1/m
    Lunit = ke0 / ke_phys                   # m per model length unit
    # amplitude: upwash variance of the e=0 spectrum = u'^2
    sp = np.load(os.path.join(HERE, "spectra_gateA_S1.npz"),
                 allow_pickle=True)

    def r22(e):
        i = int(round(e / 0.5))
        return np.trapezoid(sp[f"med_E2_e{i}"], sp[f"k2_e{i}"])

    # model-unit variance of the e=0 representation (for beta)
    kx = fl.sym_log_grid(ke0, n=fl.NKY)
    ky = fl.sym_log_grid(ke0, n=fl.NKY)
    V0m = np.trapezoid(
        [np.trapezoid([fl.phi_ww(a, b, 0.0, ke0, A00) for b in ky], ky)
         for a in kx], kx)
    beta = UP ** 2 * Lunit ** 2 / V0m       # phys amplitude per model phi

    f = np.geomspace(30.0, 2000.0, 60)
    omega = 2 * np.pi * f
    Kx = omega / U                          # physical streamwise wavenumber
    kxm = Kx * Lunit                        # model units
    kred = omega * (CHORD / 2) / U
    S2 = sears(kred)

    fig, ax = plt.subplots(figsize=(0.7 * FULL, 2.6))
    # analytic Liepmann baseline (len_vsdb form, Corcos-free: same ell_y
    # definition as the model curves would need the 2D form, so use the
    # e=0 representation's ell_y — spectrum-consistent, and at e=0 the
    # representation IS vK, so this is a pure shape cross-check)
    styles = [("frozen von Kármán (standard practice)", 0.0,
               "k", (0, (5, 2))),
              ("ODT-distorted, $e=0.5$", 0.5, COL["odt"], "-"),
              ("ODT-distorted, $e=1$", 1.0, COL["phi33"], "-"),
              ("ODT-distorted, $e=2$", 2.0, COL["rdt"], "-")]
    for lab, e, col, ls in styles:
        ke, A0 = fits[e]
        lev = r22(e) / r22(0.0)             # measured variance ratio
        p0 = np.array([fl.phi_ww(k, 0.0, e, ke, A0) for k in kxm])
        kyg = fl.sym_log_grid(ke, n=fl.NKY)
        I = np.array([np.trapezoid(
            [fl.phi_ww(k, b, e, ke, A0) for b in kyg], kyg)
            for k in kxm])
        # normalise the representation to its own variance, then apply
        # the measured level: phys Phi_ww = beta*(V0m/Ve)*lev * model
        Vem = np.trapezoid(
            [np.trapezoid([fl.phi_ww(a, b, e, ke, A0) for b in kyg],
                          kyg) for a in fl.sym_log_grid(ke, n=fl.NKY)],
            fl.sym_log_grid(ke, n=fl.NKY))
        # physical 2D density: Phi_phys(Kx,Ky) = amp * phi_model(Kx*L, Ky*L)
        # with amp = beta*(V0m/Vem)*lev; beta carries the L^2 Jacobian and
        # pins the e=0 variance to u'^2, the ratio renormalises each fitted
        # representation to unit variance, lev applies the measured level.
        # len_vsdb takes the 1-D upwash spectrum (int over Kx = w'^2):
        # Phi1D(Kx) = int Phi_phys dKy = amp*I/Lunit, and ell_y is the
        # spectrum-consistent spanwise length pi*Phi(Kx,0)/Phi1D(Kx).
        amp = beta * (V0m / Vem) * lev
        phi1d = amp * I / Lunit
        ell_y = np.pi * p0 / I * Lunit
        Spp = far_field_psd_dipole(XOBS, omega, S2, phi1d, CHORD, U,
                                   SPAN, M, ell_y)
        spl = 10 * np.log10(Spp * 2 * np.pi / PREF ** 2)   # -> dB/Hz
        ax.plot(f, spl, ls=ls, color=col, lw=1.1, label=lab)
        if e == 0.0:
            # self-validation: the e=0 representation against analytic vK
            rat = phi1d / vk_phi1d(Kx, UP ** 2, LAMBDA)
            print("e=0 Phi1D / analytic vK at 100/500/1000 Hz: "
                  + " ".join(f"{np.interp(v, f, rat):.3f}"
                             for v in (100, 500, 1000)))
            spl_lp = 10 * np.log10(far_field_psd_dipole(
                XOBS, omega, S2, liepmann_phi1d(Kx, UP ** 2, LAMBDA),
                CHORD, U, SPAN, M, ell_y) * 2 * np.pi / PREF ** 2)
            ax.plot(f, spl_lp, ls=(0, (1, 1.2)), color="0.55", lw=0.9,
                    label="frozen Liepmann")
        for ff in (100, 500, 1000):
            print(f"{lab} SPL({ff} Hz) = "
                  f"{np.interp(ff, f, spl):.1f} dB/Hz")
    ax.set_xscale("log")
    ax.set_xlabel(r"$f$ [Hz]")
    ax.set_ylabel(r"SPL [dB Hz$^{-1}$ re $20\,\mu$Pa]")
    ax.legend(fontsize=6.5, loc="lower left")
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_amiet_abs"))


if __name__ == "__main__":
    main()
