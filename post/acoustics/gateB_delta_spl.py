#!/usr/bin/env python3
r"""
Gate B: Delta-SPL of leading-edge (turbulence-interaction) noise from the
anisotropy-stretched vK upwash spectrum, relative to an isotropic vK baseline
(paper 2, Session 3).

--------------------------------------------------------------------------
Model and the two modelling calls that were fixed
--------------------------------------------------------------------------
In Amiet's leading-edge-noise theory the far-field PSD of a large-aspect-ratio
airfoil factorises as

    S_pp(x,omega) = [ rho0, M, observer geometry, semichord, |L|^2 ] (K_x)
                    * Phi_ww(K_x, 0) * ell_y(omega),

with K_x = omega/U (frozen convection), |L|^2 the airfoil gust-response
function, and the spanwise correlation length

    ell_y(omega) = pi * Phi_ww(K_x,0) / int_{-inf}^{inf} Phi_ww(K_x,K_y) dK_y.

    # CALL 1 (spanwise length): this Amiet correlation length -- the reciprocal
    #   of the note sec.5.2 ratio (axisym_family.spanwise_length) -- is the
    #   physically correct factor in S_pp and decays like 1/omega at high
    #   frequency.  The sec.5.2 ratio is U-shaped and is NOT what multiplies here.
    # CALL 2 (airfoil response): Delta-SPL compares the SAME airfoil, flow (U,M),
    #   observer, and semichord with only the turbulence spectrum changed, so all
    #   of [...] and |L|^2(K_x) cancel exactly.  No airfoil geometry or Amiet
    #   response evaluation is needed for the DELTA (it is required only for an
    #   absolute SPL, left for later).
    # Both are DECISIONS from the design chat; flip either here if note sec.5/6
    # says otherwise -- the reductions below are self-consistent regardless.

Hence the noise-relevant turbulence kernel and the level change are

    K_noise(K_x) = Phi_ww(K_x,0) * ell_y = pi * Phi_ww(K_x,0)^2 / I(K_x),
                   I(K_x) = int_{-inf}^{inf} Phi_ww(K_x,K_y) dK_y,
    Delta-SPL(omega) = 10 log10[ K_noise^aniso(K_x) / K_noise^iso(K_x) ].

Conventions and the spectral model live in axisym_family; see its docstring.
"""

import numpy as np
from scipy.special import beta

import axisym_family as af
from axisym_family import Params

try:
    _trapz = np.trapezoid
except AttributeError:  # pragma: no cover
    _trapz = np.trapz

# Log-grid for the spanwise (K_y) and K_2 integrals.  These reductions are
# smooth ratios of von Karman integrals; half=22, n=300 already matches the
# exact closed forms to ~1e-10 (see test_spanwise_corr_length_vK) at ~5x the
# speed of the axisym_family science default.
_HALF = 22.0
_N = 300


# ----------------------------------------------------------------------
# Frozen-convection map and spanwise reductions
# ----------------------------------------------------------------------
def kx_from_omega(omega, U):
    """Chordwise gust wavenumber K_x = omega / U (Taylor frozen convection)."""
    return np.asarray(omega, dtype=float) / U


def spanwise_integral(kx, p, half=_HALF, n=_N):
    """I(K_x) = int_{-inf}^{inf} Phi_ww(K_x, K_y) dK_y  (= 2 * int_0^inf)."""
    scale = max(p.L2, p.Lperp)
    ky, lnky = af._log_grid(scale, half, n)
    phiww = af.phi_ww_planar(kx, ky, p, half, n)
    return 2.0 * _trapz(phiww * ky, lnky)


def phi_ww0(kx, p, half=_HALF, n=_N):
    """Phi_ww(K_x, 0) on the same fast grid as the spanwise integral."""
    return af.phi_ww_planar(kx, 0.0, p, half, n)


def spanwise_corr_length(kx, p, half=_HALF, n=_N):
    """Amiet spanwise correlation length ell_y = pi Phi_ww(K_x,0) / I(K_x)."""
    return np.pi * phi_ww0(kx, p, half, n) / spanwise_integral(kx, p, half, n)


def noise_kernel(kx, p, half=_HALF, n=_N):
    """Turbulence factor of S_pp up to the (cancelling) airfoil/geometry terms:

        K_noise = Phi_ww(K_x,0) * ell_y = pi Phi_ww(K_x,0)^2 / I(K_x).
    """
    phi0 = phi_ww0(kx, p, half, n)
    return np.pi * phi0 * phi0 / spanwise_integral(kx, p, half, n)


# ----------------------------------------------------------------------
# Delta-SPL
# ----------------------------------------------------------------------
def delta_spl(omega, U, p_aniso, p_iso):
    """Delta-SPL(omega) = 10 log10[ K_noise^aniso / K_noise^iso ], K_x = omega/U.

    Scalar in -> float out; array omega -> array.
    """
    kx = kx_from_omega(omega, U)
    kx_arr = np.atleast_1d(kx)
    ka = np.array([noise_kernel(k, p_aniso) for k in kx_arr])
    ki = np.array([noise_kernel(k, p_iso) for k in kx_arr])
    out = 10.0 * np.log10(ka / ki)
    return float(out[0]) if np.ndim(omega) == 0 else out


def delta_spl_curve(omega, U, p_aniso, p_iso):
    """Vectorised Delta-SPL over an omega array; returns (kx, dSPL)."""
    kx = kx_from_omega(omega, U)
    return kx, delta_spl(omega, U, p_aniso, p_iso)


# ----------------------------------------------------------------------
# Isotropic baselines matched to an anisotropic family (baseline conventions)
# ----------------------------------------------------------------------
def matched_isotropic(p, mode="perp"):
    """Isotropic vK baseline for the Delta-SPL comparison.

    mode="perp" : L2 = Lperp = p.Lperp, c0 = 0, A0 unchanged (matches the
                  in-plane scale seen by the airfoil -- the natural pre-strain
                  reference for a line-normal ODT stretch).
    mode="mean" : L2 = Lperp = sqrt(p.L2*p.Lperp) (geometric-mean scale).
        # DECISION NEEDED: the pre-strain reference is a physics choice (match
        # k_t? match Lambda? match perpendicular scale?).  Exposed as `mode`;
        # confirm the intended baseline with note sec.6.
    """
    if mode == "perp":
        L = p.Lperp
    elif mode == "mean":
        L = np.sqrt(p.L2 * p.Lperp)
    else:
        raise ValueError(f"unknown mode {mode!r}")
    return Params(A0=p.A0, c0=0.0, L2=L, Lperp=L)


# ----------------------------------------------------------------------
# Isotropic-limit closed forms (for tests / eyeballing)
# ----------------------------------------------------------------------
def spanwise_corr_length_vK(kx, Lambda):
    r"""Exact isotropic vK Amiet spanwise correlation length:

        ell_y = pi k_x^2 / [ c (k_x^2 B(1/2,11/6) + B(3/2,5/6) c^2) ],
        c = sqrt(k_e^2 + k_x^2),  k_e = ke_from_length(Lambda).

    (Decays like 1/k_x at high frequency.)
    """
    ke = af.ke_from_length(Lambda)
    kx = np.asarray(kx, dtype=float)
    c = np.sqrt(ke * ke + kx * kx)
    b1 = beta(0.5, 11.0 / 6.0)
    b2 = beta(1.5, 5.0 / 6.0)
    return np.pi * kx * kx / (c * (kx * kx * b1 + b2 * c * c))


# ----------------------------------------------------------------------
# Eyeball a Delta-SPL curve when run as a script
# ----------------------------------------------------------------------
def _demo():
    U = 30.0
    Lam = 1.0
    ke = af.ke_from_length(Lam)
    # anisotropic: upwash-normal scale stretched (L2 < Lperp) -> more upwash at
    # mid/high wavenumber; keep A0 fixed and compare to the perp-matched iso.
    p_a = Params(A0=1.0, c0=0.0, L2=0.5 * ke, Lperp=ke)
    p_i = matched_isotropic(p_a, mode="perp")
    omega = np.geomspace(0.2, 40.0, 9) * U * ke  # so K_x spans ~0.2..40 * ke
    kx, d = delta_spl_curve(omega, U, p_a, p_i)
    print(f"Gate B Delta-SPL demo (U={U}, Lambda={Lam}, k_e={ke:.4f})")
    print(f"{'K_x/k_e':>9} {'omega':>10} {'dSPL[dB]':>10}")
    for kxi, oi, di in zip(kx, omega, d):
        print(f"{kxi / ke:9.3f} {oi:10.3f} {di:10.3f}")


if __name__ == "__main__":
    _demo()
