#!/usr/bin/env python3
"""SPL(f) at application strains against the DNS-anchored truth:
at the reference configuration (sec:spec-op physical bridge), rapid
rate (Sk/eps~8, where sec 4.4 shows exact RDT = DNS), compare

  - exact-RDT truth: the RDT-distorted vK of the model's own relaxed
    e=0 state, kinematic ke(e) = ke0 exp(e/2), level = the kernel's
    own exact upwash-variance amplification (checked against Level-0);
  - frozen von Karman (standard practice);
  - standard ODT and clock-relaxed ODT (fitted representations +
    measured ensemble levels), each with an SPL uncertainty band
    +-10 log10(1+rms) from its measured spectral rms distance to the
    RDT family (fig_dividend / fit costs).

Panels: e=1 and e=2.  Output fig_amiet_rapid.{pdf,png} + band table.
"""
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import rdt_kernel as rk  # noqa: E402
from figstyle_jfm import COL, GOLD, RED, FULL, panel, plt, save  # noqa: E402,E501

U, UP, LAMBDA = 20.0, 1.2, 0.055
CHORD, SPAN, XOBS = 0.386, 0.8, (0.0, 0.0, 1.0)
RHO0, C0, PREF = 1.225, 340.0, 20e-6
CVK = 0.746834
NKY, NK2 = 140, 120
# measured spectral rms distance to the RDT family (fig_dividend), S20
RMS = {"std": {1.0: 0.29, 2.0: 0.40},
       "fix": {1.0: 0.21, 2.0: 0.20}}


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


def sgrid(scale, n):
    x = scale * np.exp(np.linspace(-9.0, np.log(30.0), n // 2))
    return np.concatenate([-x[::-1], x])


def phi_m(kx, ky, e, ke):
    """model-unit Phi_ww (unnormalised A0=1)."""
    k2 = sgrid(ke, NK2)
    P22, _ = rk.rdt_phi_components(kx / ke, k2 / ke, ky / ke, e)
    return np.trapezoid(P22, k2)


def variance_m(e, ke):
    g = sgrid(ke, NKY)
    rows = np.array([np.trapezoid(
        np.array([phi_m(kx, ky, e, ke) for ky in g]), g) for kx in g])
    return np.trapezoid(rows, g)


def spl_curve(f, e, ke_m, Lunit, lev, w2):
    """SPL(f) via the Gate-B kernel chain (compact form of fig_amiet_abs:
    audited dipole + exact Sears), physical amplitude w2*lev."""
    from scipy.special import jv, yv
    omega = 2 * np.pi * f
    Kx = omega / U
    kxm = Kx * Lunit
    p0 = np.array([phi_m(k, 0.0, e, ke_m) for k in kxm])
    g = sgrid(ke_m, NKY)
    I = np.array([np.trapezoid(
        np.array([phi_m(k, ky, e, ke_m) for ky in g]), g) for k in kxm])
    V = variance_m(e, ke_m)
    amp = w2 * lev * Lunit ** 2 / V
    phi1d = amp * I / Lunit
    ell_y = np.pi * p0 / I * Lunit
    kred = omega * (CHORD / 2) / U
    H1 = jv(1, kred) - 1j * yv(1, kred)
    H0 = jv(0, kred) - 1j * yv(0, kred)
    Cth = H1 / (H1 + 1j * H0)
    S = Cth * (jv(0, kred) - 1j * jv(1, kred)) + 1j * jv(1, kred)
    M = U / C0
    sigma2 = XOBS[0] ** 2 + (1 - M ** 2) * (XOBS[1] ** 2 + XOBS[2] ** 2)
    dip2 = (omega * XOBS[2] / (4 * np.pi * C0 * sigma2)) ** 2
    lift = (np.pi * RHO0 * U * CHORD) ** 2 * np.abs(S) ** 2
    Spp = dip2 * lift * (phi1d / U) * SPAN * ell_y
    return 10 * np.log10(Spp * 2 * np.pi / PREF ** 2)


def r22(sp, e):
    i = int(round(e / 0.5))
    return np.trapezoid(sp[f"med_E2_e{i}"], sp[f"k2_e{i}"])


def main():
    fits_s = read_fits("S20")
    fits_f = read_fits("S20_RCS1")
    sp_s = np.load(os.path.join(HERE, "spectra_gateA_S20.npz"),
                   allow_pickle=True)
    sp_f = np.load(os.path.join(HERE, "spectra_gateA_S20_RCS1.npz"),
                   allow_pickle=True)
    ke0_s, ke0_f = fits_s[0.0][0], fits_f[0.0][0]
    ke_phys = CVK / LAMBDA
    w2 = UP ** 2
    f = np.geomspace(30.0, 2000.0, 40)

    # exact upwash amplification of the kernel (truth level), checked
    g22 = {e: variance_m(e, ke0_s) / variance_m(0.0, ke0_s)
           for e in (1.0, 2.0)}
    print("kernel upwash amplification G22(1), G22(2):",
          round(g22[1.0], 3), round(g22[2.0], 3),
          " (Level-0 exact at e=1: 1.498)")

    fig, axs = plt.subplots(1, 2, figsize=(FULL, 2.5), sharey=True)
    lines = []
    for ax, e in zip(axs, (1.0, 2.0)):
        # truth: kinematic ke, exact kernel level
        ke_t = ke0_s * np.exp(e / 2.0)
        spl_t = spl_curve(f, e, ke_t, ke0_s / ke_phys, g22[e], w2)
        # frozen vK
        spl_v = spl_curve(f, 0.0, ke0_s, ke0_s / ke_phys, 1.0, w2)
        # standard / fixed: fitted ke + measured level
        out = {}
        for tag, fits, sp, ke0 in (("std", fits_s, sp_s, ke0_s),
                                   ("fix", fits_f, sp_f, ke0_f)):
            lev = r22(sp, e) / r22(sp, 0.0)
            out[tag] = spl_curve(f, e, fits[e][0], ke0 / ke_phys, lev,
                                 w2)
        ax.plot(f, spl_v, "k", ls=(0, (5, 2)), lw=1.0,
                label="frozen von Kármán")
        ax.plot(f, spl_t, color="#8e44ad", lw=1.4,
                label="exact RDT (DNS-anchored)")
        for tag, col, lab in (("std", COL["odt"], "ODT standard"),
                              ("fix", RED, "ODT + clock")):
            db = 10 * np.log10(1 + RMS[tag][e])
            ax.fill_between(f, out[tag] - db, out[tag] + db, color=col,
                            alpha=0.14, lw=0)
            ax.plot(f, out[tag], color=col, lw=1.0, label=lab)
        ax.set_xscale("log")
        ax.set_xlabel(r"$f$ [Hz]")
        ax.text(0.5, 0.96, f"$e={e:.0f}$", transform=ax.transAxes,
                ha="center", va="top", fontsize=8)
        for ff in (100, 300, 1000):
            lines.append(
                f"e={e:.0f} f={ff}: truth "
                f"{np.interp(ff, f, spl_t):6.1f}  vK "
                f"{np.interp(ff, f, spl_v):6.1f}  std "
                f"{np.interp(ff, f, out['std']):6.1f}  fix "
                f"{np.interp(ff, f, out['fix']):6.1f}")
    axs[0].set_ylabel(r"SPL [dB Hz$^{-1}$ re $20\,\mu$Pa]")
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, ncol=4, loc="lower center", fontsize=6.5,
               bbox_to_anchor=(0.5, -0.12))
    panel(axs[0], "a")
    panel(axs[1], "b")
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_amiet_rapid"))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
