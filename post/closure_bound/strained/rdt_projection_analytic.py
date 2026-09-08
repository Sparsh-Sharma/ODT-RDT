#!/usr/bin/env python3
"""Analytic exact-RDT line projection: transverse splitting phi11/phi33 - 1
on the x2 line from the closed-form (Cauchy) plane-strain RDT of
lenoise/rdt_kernel.py, with an isotropic von Karman initial spectrum.

Purpose: replace the realisation-based "exact RDT" curve of
rdt_projection.npz (4 seeds of the 128^3 pilot; one integer kappa_2 line per
low-wavenumber band, hence jagged) by the deterministic integral it estimates.

    phi_ii(k2; e) = \\iint Phi_ii(k1, k2, k3; e) dk1 dk3
                  = \\int d ln kp  2 pi kp^2  < Phi_ii >_phi(kp, k2)

Wavenumber bridge to the DNS: the vK scale ke is fixed so that the e = 0
total line spectrum has centroid KC = 3 (box units) over the DNS band
(0, KMAX], exactly the normalisation of rdt_projection.py (KREF = 3).

Writes rdt_projection_analytic.npz (x = k2(e)/k_c(0), split_e0.5, split_e1,
R_nn as well for completeness) and an overlay PNG against the 4-seed curve.
"""
import os
import sys

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
LEN = os.path.normpath(os.path.join(HERE, "..", "lenoise"))
sys.path.insert(0, LEN)
import axisym_family as af  # noqa: E402
import rdt_kernel as rk  # noqa: E402

NKP, NPHI = 600, 128
KC = 3.0              # e=0 centroid of the total line spectrum (box units)
KMAX = 0.85 * 42.0    # DNS dealiasing cutoff on kappa_2 (box units)


def phi_diag(k1, k2, k3, e, band=None):
    """Exact plane-strain RDT diagonal Phi_11, Phi_22, Phi_33 (same algebra
    as rdt_kernel.rdt_phi_components, all three components returned).

    band=(lo, hi): keep only initial wavevectors lo <= |k0| <= hi (the
    discrete box has no mode below the first integer wavenumber and none
    above its dealiasing cutoff) — the DNS-companion normalisation."""
    beta = np.exp(0.5 * e)
    f = np.array([beta, 1.0 / beta, 1.0])
    k1, k2, k3 = np.broadcast_arrays(np.asarray(k1, float),
                                     np.asarray(k2, float),
                                     np.asarray(k3, float))
    k = np.stack([k1, k2, k3])
    ksq = (k ** 2).sum(axis=0)
    k0 = k * f.reshape((3,) + (1,) * (k.ndim - 1))
    k0mag = np.sqrt((k0 ** 2).sum(axis=0))
    q = rk._E0(k0mag) / (4.0 * np.pi)
    if band is not None:
        q = q * ((k0mag >= band[0]) & (k0mag <= band[1]))
    eps = np.zeros((3, 3, 3))
    eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1.0
    eps[0, 2, 1] = eps[2, 1, 0] = eps[1, 0, 2] = -1.0
    M = np.einsum("ipa,p...,a->ia...", eps, k, f) / ksq
    k0hat = k0 / k0mag
    MMT = np.einsum("ia...,ja...->ij...", M, M)
    Mk = np.einsum("ia...,a...->i...", M, k0hat)
    Phi = q * (MMT - np.einsum("i...,j...->ij...", Mk, Mk))
    return Phi[0, 0], Phi[1, 1], Phi[2, 2]


def line_spectra(k2hat, e, nkp=NKP, nphi=NPHI, band=None):
    """phi_11, phi_22, phi_33 on the x2 line, normalised units (ke = 1)."""
    kp, lnkp = af._log_grid(1.0, 14.0, nkp)
    phi = (np.arange(nphi) + 0.5) * (2.0 * np.pi / nphi)
    K1 = kp[None, :, None] * np.cos(phi)[None, None, :]
    K3 = kp[None, :, None] * np.sin(phi)[None, None, :]
    K2 = np.asarray(k2hat, float)[:, None, None] * np.ones_like(K1)
    P11, P22, P33 = phi_diag(K1, K2, K3, e, band=band)
    w = 2.0 * np.pi * kp * kp
    out = [np.trapezoid(P.mean(axis=2) * w[None, :], lnkp, axis=1)
           for P in (P11, P22, P33)]
    return out


KLO = 1.0             # first box mode (box units)


def box_band(ke):
    """The DNS box's initial-wavevector band in normalised (ke = 1) units."""
    return (KLO / ke, KMAX / ke)


def fit_ke():
    """ke (box units) such that the e=0 total line spectrum, restricted to
    the box band, has centroid KC over (0, KMAX] — the DNS normalisation."""
    def centroid_minus_kc(ke):
        # coarse grids: the centroid is a smooth scalar, and at e=0 the
        # field is isotropic so the azimuthal average is exact at any nphi
        k2 = np.geomspace(1e-3, KMAX / ke, 120)          # normalised
        E = sum(line_spectra(k2, 0.0, nkp=160, nphi=8, band=box_band(ke)))
        c_hat = np.trapezoid(k2 * E, k2) / np.trapezoid(E, k2)
        return ke * c_hat - KC
    return brentq(centroid_minus_kc, 0.2, 6.0, xtol=1e-4)


def main():
    ke = fit_ke()
    print(f"ke (box units) matching e=0 centroid {KC} in the box band "
          f"[{KLO}, {KMAX:.1f}]: {ke:.3f}")
    x = np.geomspace(0.4, 12.0, 60)                    # k2(e)/k_c(0)
    k2hat = x * KC / ke
    res = {"x": x, "ke": ke}
    band = box_band(ke)
    # e=0 normalisation over the resolved band (0, KMAX] (box units)
    kg0 = np.geomspace(1e-2, KMAX / ke, 150)
    N0 = [np.trapezoid(p, kg0) for p in line_spectra(kg0, 0.0, band=band)]
    for e, tag in ((0.5, "0.5"), (1.0, "1")):
        f = np.exp(0.5 * e)
        P = line_spectra(k2hat, e, band=band)
        res[f"split_e{tag}"] = P[0] / P[2] - 1.0
        # shape ratio R_nn against the rigid translation k2 -> k2 f
        # (same construction as rdt_projection.py: normalised strained
        # spectrum over normalised initial spectrum at k2/f, times 1/f)
        P0 = line_spectra(k2hat / f, 0.0, band=band)
        kge = np.geomspace(1e-2, KMAX * f / ke, 150)
        Ne = [np.trapezoid(p, kge) for p in line_spectra(kge, e, band=band)]
        for c, name in ((0, "phi11"), (1, "phi22"), (2, "phi33")):
            s1 = P[c] / Ne[c]
            s0 = P0[c] / N0[c]
            res[f"R_{name}_e{tag}"] = s1 / (s0 / f)
    np.savez(os.path.join(HERE, "rdt_projection_analytic.npz"), **res)

    # ---- overlay against the 4-seed companions ------------------------
    D = np.load(os.path.join(HERE, "rdt_projection.npz"), allow_pickle=True)
    XC = np.sqrt(D["edges"][:-1] * D["edges"][1:])
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5, 3.4))
    for tag, lw in (("1", 1.6), ("0.5", 0.9)):
        ax.plot(x, res[f"split_e{tag}"], "k-", lw=lw,
                label=f"analytic exact RDT, e={tag}")
        ax.plot(XC, D[f"split_rdt_e{tag}"], "o--", color="0.5", lw=lw,
                ms=4, mfc="none", label=f"4-seed 128^3 companions, e={tag}")
        # numeric comparison at the band centres
        interp = np.interp(XC, x, res[f"split_e{tag}"])
        print(f"\ne={tag}:  x      analytic   4-seed    diff")
        for xc, a_, s_ in zip(XC, interp, D[f"split_rdt_e{tag}"]):
            print(f"      {xc:5.2f}   {a_:+.3f}    {s_:+.3f}   {a_ - s_:+.3f}")
        for name in ("phi22", "phi11", "phi33"):
            ia = np.interp(XC, x, res[f"R_{name}_e{tag}"])
            print(f"  R_{name}: analytic "
                  + " ".join(f"{v:5.2f}" for v in ia))
            print(f"  R_{name}: 4-seed   "
                  + " ".join(f"{v:5.2f}" for v in D[f"R_{name}_e{tag}"]))
    ax.axhline(0, color="0.8", lw=0.6, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("k2(e)/k_c(0)")
    ax.set_ylabel("phi11/phi33 - 1")
    ax.legend(fontsize=7)
    fig.tight_layout()
    out = os.path.join(HERE, "rdt_projection_analytic_check.png")
    fig.savefig(out, dpi=150)
    print("saved", out)


if __name__ == "__main__":
    main()
