#!/usr/bin/env python3
r"""Angular kernel h(mu; e) for the axisym family, derived from exact
plane-strain rapid distortion theory (replaces the P2-Legendre stub, which
the gateA_S1 fits showed to be structurally insufficient at e >= 1).

Derivation
----------
Plane strain A = diag(a, -a, 0), total strain e with the Test-3 convention
e = S t, S = 1, a = 1/2, so the principal material stretch is
beta = exp(a t) = exp(e/2):  F = diag(beta, 1/beta, 1).

Exact (Cauchy) RDT per mode, starting from isotropic turbulence:

    k0    = F^T k                     (wavevector back-map)
    omega = F omega0                  (Cauchy vorticity)
    u_i   = i eps_ipq k_p omega_q / k^2

    <omega0_a omega0_b*>(k0) = Omega0_ab = (E0(k0)/4pi) (delta_ab - k0hat_a k0hat_b)

    Phi_ij(k; e) = M_ia M_jb Omega0_ab,   M_ia = eps_ipa k_p f_a / k^2
    (F = diag(f1,f2,f3); no sum over the diagonal index a in f_a).

This Phi is NOT axisymmetric about e_2 (plane strain has three axes); the
axisymmetric family can only represent its RING AVERAGE about e_2.  With the
family's ring reductions (axisym_family._Phi_ii_ring)

    4pi kap^2 <Phi_22>_ring          = A (1-mu^2) + C (1-mu^2)^2
    4pi kap^2 <Phi_11 + Phi_33>_ring = A (1+mu^2) + C mu^2 (1-mu^2)

the pair (A_eff, C_eff) is an exact 2x2 inversion at every (mu, kap), and

    h_RDT(mu; e) = C_eff / A_eff   (normalized to max|h| = 1; the family's
                                    c0 carries the amplitude).

In the inertial (power-law) band of E0 the RDT map is scale-free, so
h_RDT is kappa-independent there — verified by the collapse check in the
demo.  E0 is the vK shape from axisym_family (ke = 1).

Usage: h = kernel_for_e(e) gives a callable h(mu) for axisym_family's
set_h_kernel; `python3 rdt_kernel.py` runs the isotropy test, the collapse
check, and writes fig_rdt_kernel.png + rdt_kernel_tab.npz.
"""
import os

import numpy as np

import axisym_family as af

NPHI = 256          # azimuth nodes for the ring average
KAPPA_REF = 30.0    # inertial-band reference kappa/ke for the kernel table
MU = np.linspace(0.0, 0.999, 121)


def _E0(k):
    """Isotropic vK energy spectrum shape, ke = 1 (amplitude irrelevant)."""
    return af.Psi(k)


def rdt_phi_components(k1, k2, k3, e):
    """Exact plane-strain RDT Phi_22 and Phi_11+Phi_33 at wavevector(s) k.

    Broadcasts over arrays.  Amplitude carries E0's normalization.
    """
    beta = np.exp(0.5 * e)
    f = np.array([beta, 1.0 / beta, 1.0])
    k1, k2, k3 = np.broadcast_arrays(np.asarray(k1, float),
                                     np.asarray(k2, float),
                                     np.asarray(k3, float))
    k = np.stack([k1, k2, k3])                       # (3, ...)
    ksq = (k ** 2).sum(axis=0)
    k0 = k * f[:, None] if k.ndim == 2 else k * f.reshape((3,) + (1,) * (k.ndim - 1))
    k0sq = (k0 ** 2).sum(axis=0)
    k0mag = np.sqrt(k0sq)
    q = _E0(k0mag) / (4.0 * np.pi)                   # Omega0 scalar

    # M_ia = eps_ipa k_p f_a / k^2   (i rows, a cols)
    eps = np.zeros((3, 3, 3))
    eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1.0
    eps[0, 2, 1] = eps[2, 1, 0] = eps[1, 0, 2] = -1.0
    M = np.einsum("ipa,p...,a->ia...", eps, k, f) / ksq

    # Phi_ij = M_ia M_jb q (delta_ab - k0hat_a k0hat_b)
    k0hat = k0 / k0mag
    MMT = np.einsum("ia...,ja...->ij...", M, M)
    Mk = np.einsum("ia...,a...->i...", M, k0hat)
    Phi = q * (MMT - np.einsum("i...,j...->ij...", Mk, Mk))
    return Phi[1, 1], Phi[0, 0] + Phi[2, 2]


def ring_averages(mu, kap, e, nphi=NPHI):
    """Ring-averaged (Phi22, Phi11+Phi33) about e_2 at fixed (mu, kap)."""
    mu = np.atleast_1d(np.asarray(mu, float))
    kperp = kap * np.sqrt(np.clip(1.0 - mu * mu, 0.0, None))
    phi = (np.arange(nphi) + 0.5) * (2.0 * np.pi / nphi)
    k1 = kperp[:, None] * np.cos(phi)[None, :]
    k3 = kperp[:, None] * np.sin(phi)[None, :]
    k2 = (kap * mu)[:, None] * np.ones_like(k1)
    P22, Pperp = rdt_phi_components(k1, k2, k3, e)
    return P22.mean(axis=1), Pperp.mean(axis=1)


def invert_AC(mu, kap, e, nphi=NPHI):
    """Exact (A_eff, C_eff) of the axisym family reproducing the RDT rings."""
    r22, rpp = ring_averages(mu, kap, e, nphi)
    s = 4.0 * np.pi * kap * kap
    y1, y2 = s * r22, s * rpp
    om = 1.0 - mu * mu
    # y1 = A om + C om^2 ; y2 = A (1+mu^2) + C mu^2 om
    det = om * mu * mu * om - om * om * (1.0 + mu * mu)
    A = (y1 * mu * mu * om - om * om * y2) / det
    C = (om * y2 - (1.0 + mu * mu) * y1) / det
    return A, C


def kernel_table(e, kap=KAPPA_REF, mu=MU):
    """h_RDT(mu; e) normalized to max|h| = 1, plus the amplitude c0_eff
    (median of C/A over mu, against the normalized shape)."""
    A, C = invert_AC(mu, kap, e)
    h = C / A
    scale = np.max(np.abs(h))
    if scale == 0.0:
        return np.zeros_like(mu), 0.0
    hn = h / scale
    return hn, scale


def kernel_for_e(e, kap=KAPPA_REF):
    """Callable h(mu) (even in mu, clipped ends) for af.set_h_kernel."""
    hn, _ = kernel_table(e, kap)
    def h(mu):
        m = np.clip(np.abs(np.asarray(mu, float)), MU[0], MU[-1])
        return np.interp(m, MU, hn)
    return h


# ----------------------------------------------------------------------
def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # isotropy limit: C/A -> 0 at e = 0
    A0_, C0_ = invert_AC(MU[1:-1], KAPPA_REF, 0.0)
    print(f"e=0 isotropy check: max|C/A| = {np.max(np.abs(C0_/A0_)):.2e}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.2))
    for e, col in zip((0.5, 1.0, 1.5, 2.0), ("C0", "C1", "C3", "C2")):
        for kap, ls in ((10.0, ":"), (30.0, "-"), (100.0, "--")):
            hn, sc = kernel_table(e, kap)
            ax1.plot(MU, hn, ls, color=col, lw=1.1,
                     label=f"e={e}, kap={kap:.0f}" if kap == 30 else None)
        A, C = invert_AC(MU, KAPPA_REF, e)
        ax2.plot(MU, C / A, color=col, label=f"e={e}: c0_eff={sc:.2f}")
    ax1.plot(MU, 0.5 * (3 * MU ** 2 - 1), "k-.", lw=1.4, label="P2 stub")
    ax1.set_xlabel(r"$\mu$"); ax1.set_ylabel(r"$h(\mu)$ (normalized)")
    ax1.set_title("RDT kernel: shape and $\\kappa$-collapse")
    ax1.legend(fontsize=7)
    ax2.set_xlabel(r"$\mu$"); ax2.set_ylabel(r"$C_{eff}/A_{eff}$")
    ax2.set_title("unnormalized (amplitude = family's $c_0$)")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    here = os.path.dirname(os.path.abspath(__file__))
    fig.savefig(os.path.join(here, "results_test3", "fig_rdt_kernel.png"), dpi=180)
    tabs = {f"h_e{e}": kernel_table(e)[0] for e in (0.5, 1.0, 1.5, 2.0)}
    np.savez(os.path.join(here, "results_test3", "rdt_kernel_tab.npz"),
             mu=MU, **tabs)
    print("saved fig_rdt_kernel.png + rdt_kernel_tab.npz")


if __name__ == "__main__":
    _demo()
