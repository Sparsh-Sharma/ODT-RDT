#!/usr/bin/env python3
r"""
Parametric axisymmetric upwash-spectrum family (anisotropy-stretched von Karman)
and its projections, for paper 2 (Gate A).

This module implements the four-parameter solenoidal axisymmetric velocity
spectrum family and the reductions needed for the upwash (w = u2) statistics
that feed an Amiet-style airfoil-turbulence-interaction response. It contains
NO fitting and reads NO run data -- it is pure spectral math plus the closed
forms it must reduce to.  See paper2_gateA_derivation.md sec.3, sec.5, sec.6.

--------------------------------------------------------------------------
Conventions (must match paper 1)
--------------------------------------------------------------------------
Axes            : x1 chordwise, x2 airfoil-normal (= ODT line y), x3 spanwise.
Wavevector      : kappa = (k1, k2, k3);  k2 is the ODT-resolved axis.
                  kperp = sqrt(k1^2 + k3^2);  kap = sqrt(k1^2+k2^2+k3^2).
Symmetry axis   : e = e_2 (the ODT/airfoil-normal direction).  The field is
                  axisymmetric about e_2, so mu = (e . kappa)/kap = k2/kap.
Amiet maps      : k_x = k1 (chordwise), k_y = k3 (spanwise).
Upwash          : w = u2; its spectral density is Phi_22.
Energy norm     : int E_i(k2) dk2 = R_ii (component variance, no sum), and
                  (1/2) sum_i R_ii = k_t.

--------------------------------------------------------------------------
The family (note sec.3.1)
--------------------------------------------------------------------------
Anisotropy-stretched von Karman shape over (k2, kperp):

    khat^2   = (k2/L2)^2 + (kperp/Lperp)^2
    Psi(khat)= khat^4 / (1 + khat^2)^(17/6)          # von Karman shape
    A(k2,kp) = A0 * Psi(khat)                         # ~ energy spectrum E(k)
    C(k2,kp) = c0 * A(k2,kp) * h(mu),  mu = k2/kap

The solenoidal axisymmetric tensor is built by incompressible projection of the
axisymmetric "pre-tensor"  T_kl = A*delta_kl + C*e_k e_l  (e = e_2):

    Phi_ij(kappa) = (1/(4 pi kap^2)) * P_ik P_jl T_kl
                  = (1/(4 pi kap^2)) * ( A * P_ij  +  C * eperp_i eperp_j )

with the isotropic projector  P_ij = delta_ij - k_i k_j / kap^2  and the
projected symmetry axis  eperp_i = P_ij e_j = e_i - mu * k_i/kap.
This is manifestly solenoidal (k_i Phi_ij = 0) and, for c0 = 0 and
L2 = Lperp, collapses to the isotropic von Karman tensor
Phi_ij = E(k)/(4 pi k^2) (delta_ij - k_i k_j/k^2) with E(k) = A0*Psi.

    # DECISION NEEDED: the second solenoidal structure is taken here as the
    # rank-1 projected-axis dyad eperp_i eperp_j (i.e. the incompressible
    # projection of C * e_i e_j).  paper2_gateA_derivation.md sec.3 is not
    # reachable in this repo; if note sec.3 uses a different (e.g. traceless
    # Q_ij = eperp_i eperp_j - (1/2)|eperp|^2 P_ij) C-tensor, only the c0 != 0
    # branch changes -- the isotropic reductions (tests 1,2,3,5) are unaffected
    # because C == 0 there.  Paste this back to the design chat to confirm.

--------------------------------------------------------------------------
Length- vs wavenumber-scale convention for L2, Lperp
--------------------------------------------------------------------------
khat = kappa / L, so L2, Lperp are ENERGY-WAVENUMBER scales (units of a
wavenumber), NOT lengths -- this keeps khat dimensionless and makes Psi the
standard (k/k_e) von Karman shape.  For an isotropic field of integral length
scale Lambda set L2 = Lperp = ke_from_length(Lambda), with

    k_e = sqrt(pi) * Gamma(5/6) / (Lambda * Gamma(1/3)).

The isotropic-limit closed forms in sec.6 use exactly this k_e, so the shape
test (phi_ww_ky0) matches only under this reading of L; see ke_from_length.

--------------------------------------------------------------------------
Numerics
--------------------------------------------------------------------------
The 3-D energy integrals (E_component, R_ii) are done on uniform grids in
log-wavenumber: the von Karman integrands (~ kappa^{-11/3} tail, finite core)
become smooth, exponentially decaying bumps in ln(kappa), for which the
trapezoidal rule is spectrally accurate (Euler-Maclaurin).  Default grid
LOG_HALF_WIDTH = 30 decades-ish half-width, N_LOG = 800 points, gives the
component energies to ~1e-9 (see test_parseval_R22).  The 1-D reductions
(phi_ww_planar, spanwise_length) use adaptive scipy.integrate.quad.
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.special import gamma, beta

# np.trapezoid (NumPy >= 2.0) with np.trapz fallback (match paper1 habit).
try:
    _trapz = np.trapezoid
except AttributeError:  # pragma: no cover - older numpy
    _trapz = np.trapz

# Log-grid defaults for the 3-D (E_component / R_ii) integrals.
LOG_HALF_WIDTH = 30.0  # half-width in ln(kappa) about the energy scale
N_LOG = 800            # number of log-spaced nodes


# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
@dataclass
class Params:
    """The four family parameters.

    A0    : overall amplitude of the energy spectrum (sets u'^2 / k_t).
    c0    : anisotropy weight of the C (projected-axis) contribution; c0 = 0
            gives the isotropic-in-C limit.
    L2    : energy-wavenumber scale along the ODT axis k2 (see module docstring).
    Lperp : energy-wavenumber scale in the perpendicular (k1,k3) plane.
    """

    A0: float = 1.0
    c0: float = 0.0
    L2: float = 1.0
    Lperp: float = 1.0


def ke_from_length(Lambda):
    """von Karman energy wavenumber k_e from the integral length scale Lambda.

        k_e = sqrt(pi) * Gamma(5/6) / (Lambda * Gamma(1/3)).
    """
    return np.sqrt(np.pi) * gamma(5.0 / 6.0) / (Lambda * gamma(1.0 / 3.0))


def isotropic_params(Lambda, A0=1.0):
    """Params for an isotropic von Karman field of integral scale Lambda."""
    ke = ke_from_length(Lambda)
    return Params(A0=A0, c0=0.0, L2=ke, Lperp=ke)


# ----------------------------------------------------------------------
# Scalar building blocks
# ----------------------------------------------------------------------
def Psi(khat):
    """von Karman shape  Psi = khat^4 / (1 + khat^2)^(17/6)."""
    kh2 = np.asarray(khat, dtype=float) ** 2
    return kh2 * kh2 / np.power(1.0 + kh2, 17.0 / 6.0)


def _khat2(k2, kperp, p):
    return (k2 / p.L2) ** 2 + (kperp / p.Lperp) ** 2


def A_scalar(k2, kperp, p):
    """Isotropic-like scalar A(k2, kperp) = A0 * Psi(khat)."""
    return p.A0 * Psi(np.sqrt(_khat2(k2, kperp, p)))


def _h_mu(mu):
    """Lowest-order angular kernel h(mu).

    # CONVENTION: matches paper1 g-kernels (rapid pressure-strain g1,g2,g3).
    # paper1 sec.5 kernel source (section5_rapid_pressure_strain.tex / kernel
    # code) is NOT present in this repo, so this is the P2-Legendre stub
    # h(mu) = (3 mu^2 - 1)/2.
    # TODO confirm against paper1 (cite the g-kernel equation number once the
    # source is reachable); do not silently swap in a different form.
    """
    mu = np.asarray(mu, dtype=float)
    return 0.5 * (3.0 * mu * mu - 1.0)


# Active angular kernel.  Default = the P2 stub above; swap in a derived
# kernel (e.g. rdt_kernel.kernel_for_e(e)) with set_h_kernel and restore
# with set_h_kernel(None).  The kernel must be even in mu and vectorized.
H_KERNEL = _h_mu


def set_h_kernel(h=None):
    """Install h(mu) as the family's angular kernel (None -> P2 stub)."""
    global H_KERNEL
    H_KERNEL = _h_mu if h is None else h


def C_scalar(k2, kperp, p):
    """Anisotropy scalar C(k2, kperp) = c0 * A(k2, kperp) * h(mu)."""
    k2 = np.asarray(k2, dtype=float)
    kperp = np.asarray(kperp, dtype=float)
    kap = np.hypot(k2, kperp)
    mu = np.divide(k2, kap, out=np.zeros_like(kap), where=kap > 0.0)
    return p.c0 * A_scalar(k2, kperp, p) * H_KERNEL(mu)


# ----------------------------------------------------------------------
# Spectral tensor: upwash density Phi_22, and the diagonal Phi_ii
# ----------------------------------------------------------------------
def Phi22(k1, k2, k3, p):
    """Upwash spectral density Phi_22 at a kappa point (broadcasts over arrays).

        Phi_22 = (1/(4 pi kap^2)) ( A*(1-mu^2) + C*(1-mu^2)^2 ),
        mu = k2/kap,  1 - mu^2 = kperp^2 / kap^2.
    """
    k1 = np.asarray(k1, dtype=float)
    k2 = np.asarray(k2, dtype=float)
    k3 = np.asarray(k3, dtype=float)
    kperp2 = k1 * k1 + k3 * k3
    kap2 = kperp2 + k2 * k2
    safe = kap2 > 0.0
    inv_kap2 = np.divide(1.0, kap2, out=np.zeros_like(kap2 * 1.0), where=safe)
    one_m_mu2 = kperp2 * inv_kap2  # = 1 - (k2/kap)^2 = eperp_2
    kperp = np.sqrt(kperp2)
    A = A_scalar(k2, kperp, p)
    C = C_scalar(k2, kperp, p)
    val = (A * one_m_mu2 + C * one_m_mu2 * one_m_mu2) * inv_kap2 / (4.0 * np.pi)
    if np.ndim(val) == 0:
        return float(val)
    return val


def _Phi_ii_ring(i, k2, kperp, p):
    """Azimuth-averaged Phi_ii at fixed (k2, kperp), averaged over the azimuth
    of (k1, k3).  Uses <k1^2>_ring = <k3^2>_ring = kperp^2/2 and the field's
    axisymmetry about e_2.  i in {1,2,3}; broadcasts over array (k2, kperp).

        Phi_ii = (1/(4 pi kap^2)) ( A*<P_ii> + C*<eperp_i^2> ).
    """
    k2 = np.asarray(k2, dtype=float)
    kperp = np.asarray(kperp, dtype=float)
    kap2 = k2 * k2 + kperp * kperp
    safe = kap2 > 0.0
    inv_kap2 = np.divide(1.0, kap2, out=np.zeros_like(kap2), where=safe)
    mu2 = (k2 * k2) * inv_kap2
    one_m_mu2 = (kperp * kperp) * inv_kap2
    A = A_scalar(k2, kperp, p)
    C = C_scalar(k2, kperp, p)
    if i == 2:
        # e_2 = 1: P_22 = 1 - mu^2 ; eperp_2 = 1 - mu^2.
        Pii = one_m_mu2
        eperp2 = one_m_mu2 * one_m_mu2
    else:
        # i in {1,3}: e_i = 0.  <k_i^2>_ring = kperp^2/2.
        #   <P_ii>     = 1 - <k_i^2>/kap^2 = 1 - kperp^2/(2 kap^2)
        #   eperp_i    = -mu k_i/kap  ->  <eperp_i^2> = mu^2 <k_i^2>/kap^2
        half_one_m_mu2 = 0.5 * one_m_mu2
        Pii = 1.0 - half_one_m_mu2
        eperp2 = mu2 * half_one_m_mu2
    return (A * Pii + C * eperp2) * inv_kap2 / (4.0 * np.pi)


# ----------------------------------------------------------------------
# Log-wavenumber grids for the 3-D integrals
# ----------------------------------------------------------------------
def _log_grid(scale, half=LOG_HALF_WIDTH, n=N_LOG):
    """Uniform grid in ln(kappa) centred on `scale`, returned as (kappa, lnk)."""
    lnk = np.linspace(np.log(scale) - half, np.log(scale) + half, n)
    return np.exp(lnk), lnk


def default_k2_grid(p, half=LOG_HALF_WIDTH, n=N_LOG):
    """Log-spaced positive k2 grid suitable for E_component / R_ii.

    Uniform in ln(k2) about the energy scale so the trapezoidal integrals in
    R_ii are spectrally accurate.
    """
    scale = max(p.L2, p.Lperp)
    return _log_grid(scale, half, n)[0]


# ----------------------------------------------------------------------
# One-dimensional line spectra E_i(k2) and component energies R_ii
# ----------------------------------------------------------------------
def E_component(i, k2_grid, p, half=LOG_HALF_WIDTH, n=N_LOG):
    r"""1-D line spectrum E_i(k2) for i in {1,2,3} over a grid of k2 values.

        E_i(k2) = int int Phi_ii(k1,k2,k3) dk1 dk3        (no sum on i)
                = int_0^inf 2 pi kperp <Phi_ii>_ring dkperp

    By axisymmetry about e_2 the (k1,k3) integral collapses to the radial
    integral above; it is evaluated on a log-kperp grid (spectrally accurate
    for the von Karman integrand).  Returns an array shaped like k2_grid
    (scalar in -> float out).
    """
    scale = max(p.L2, p.Lperp)
    kperp, lnkp = _log_grid(scale, half, n)          # (n,)
    k2_arr = np.atleast_1d(np.asarray(k2_grid, dtype=float))  # (m,)
    K2 = k2_arr[:, None]                              # (m,1)
    KP = kperp[None, :]                               # (1,n)
    ring = _Phi_ii_ring(i, K2, KP, p)                # (m,n)
    # int f dkperp = int (f*kperp) dln(kperp)
    integrand = 2.0 * np.pi * KP * ring * KP         # 2 pi kperp * ring * kperp
    E = _trapz(integrand, lnkp, axis=1)              # (m,)
    if np.ndim(k2_grid) == 0:
        return float(E[0])
    return E


def R_ii(k2_grid, p):
    """Component energies [R_11, R_22, R_33] = int E_i(k2) dk2  (no sum on i).

    k2_grid is a POSITIVE, log-spaced half-axis (see default_k2_grid); each
    E_i is even in k2 so R_ii = 2 * int_0^inf E_i dk2, evaluated as
    2 * int E_i(k2) k2 dln(k2) by the trapezoidal rule.  (1/2) sum_i R_ii = k_t.
    """
    k2_grid = np.asarray(k2_grid, dtype=float)
    lnk2 = np.log(k2_grid)
    out = np.empty(3)
    for idx, i in enumerate((1, 2, 3)):
        Ei = E_component(i, k2_grid, p)
        out[idx] = 2.0 * _trapz(Ei * k2_grid, lnk2)
    return out


# ----------------------------------------------------------------------
# Amiet-plane upwash spectrum and spanwise correlation length (note sec.5)
# ----------------------------------------------------------------------
def phi_ww_planar(kx, ky, p, half=LOG_HALF_WIDTH, n=N_LOG):
    r"""Planar upwash spectrum  Phi_ww(k_x,k_y) = int Phi_22(kx,k2,ky) dk2.

    k_x = k1 (chordwise), k_y = k3 (spanwise); k2 (airfoil-normal) is integrated
    out.  Phi_22 is even in k2, so the integral is 2 * int_0^inf, evaluated on a
    log-k2 grid.  Broadcasts over array-valued ky (kx scalar).
    """
    scale = max(p.L2, p.Lperp)
    k2, lnk2 = _log_grid(scale, half, n)             # (n,)
    ky_arr = np.atleast_1d(np.asarray(ky, dtype=float))
    P = Phi22(kx, k2[None, :], ky_arr[:, None], p)   # (m,n)
    out = 2.0 * _trapz(P * k2[None, :], lnk2, axis=1)  # int = int Phi*k2 dln k2
    if np.ndim(ky) == 0:
        return float(out[0])
    return out


def phi_ww_ky0(kx, p):
    """Planar upwash spectrum on the k_y = 0 line, Phi_ww(k_x, 0)."""
    return phi_ww_planar(kx, 0.0, p)


def spanwise_length(kx, p, half=LOG_HALF_WIDTH, n=N_LOG):
    r"""Spanwise length from the note sec.5.2 definition

        ell_y(k_x) = (1/Phi_ww(k_x,0)) int_0^inf Phi_ww(k_x,k_y) dk_y.

    NOTE (non-monotone): with this literal definition the isotropic von Karman
    closed form is U-shaped in k_x (see spanwise_length_vK), diverging as
    k_x -> 0 and k_x -> inf with a single interior minimum -- it is NOT monotone
    decreasing.  The classical Amiet spanwise *correlation* length is the
    reciprocal-like quantity  pi*Phi_ww(k_x,0)/int Phi_ww dk_y, which does
    decrease with k_x.
        # DECISION NEEDED: the session brief expects a monotone-decreasing
        # trend, which contradicts the closed form of this exact ratio.  Confirm
        # with note sec.5.2 whether ell_y is this ratio (verified below vs the
        # exact closed form) or the reciprocal Amiet correlation length.
    """
    denom = phi_ww_ky0(kx, p)
    ky, lnky = _log_grid(max(p.L2, p.Lperp), half, n)  # (m,)
    phiww = phi_ww_planar(kx, ky, p, half, n)          # (m,)
    # int_0^inf Phi_ww dky = int Phi_ww * ky dln(ky)
    num = _trapz(phiww * ky, lnky)
    return num / denom


# ----------------------------------------------------------------------
# Isotropic-limit closed forms (note sec.6) -- for tests / eyeballing
# ----------------------------------------------------------------------
def phi_ww_ky0_vK_shape(kx, Lambda):
    """von Karman transverse upwash spectrum shape on k_y = 0 (up to a constant):

        Phi_ww^vK(k_x,0) ~ (k_x/k_e)^2 / (1 + (k_x/k_e)^2)^(7/3),
        k_e = sqrt(pi) Gamma(5/6)/(Lambda Gamma(1/3)).
    """
    ke = ke_from_length(Lambda)
    r2 = (np.asarray(kx, dtype=float) / ke) ** 2
    return r2 / np.power(1.0 + r2, 7.0 / 3.0)


def spanwise_length_vK(kx, Lambda):
    r"""Exact isotropic von Karman closed form of the sec.5.2 spanwise length

        ell_y(k_x) = (1/Phi_ww(k_x,0)) int_0^inf Phi_ww(k_x,k_y) dk_y

    with the 2-D upwash spectrum Phi_ww(k_x,k_y) ~ (k_x^2+k_y^2)/
    (k_e^2+k_x^2+k_y^2)^(7/3).  Carrying out the k_y integral (Beta integrals):

        ell_y = (1/2) B(1/2,11/6) * c
              + (1/2) B(3/2,5/6) * c^3 / k_x^2,   c = sqrt(k_e^2 + k_x^2).

    (U-shaped: -> inf as k_x -> 0 and k_x -> inf.)
    """
    ke = ke_from_length(Lambda)
    kx = np.asarray(kx, dtype=float)
    c = np.sqrt(ke * ke + kx * kx)
    return (0.5 * beta(0.5, 11.0 / 6.0) * c
            + 0.5 * beta(1.5, 5.0 / 6.0) * c ** 3 / (kx * kx))


# ----------------------------------------------------------------------
# Eyeball overlay when run as a script (Definition of done, sec. print block)
# ----------------------------------------------------------------------
def _print_overlay(Lambda=1.0, A0=1.0):
    p = isotropic_params(Lambda, A0=A0)
    ke = ke_from_length(Lambda)
    kxs = np.array([0.25, 0.5, 1.0, 2.0, 4.0]) * ke
    print(f"\nIsotropic reduction overlay (Lambda={Lambda}, A0={A0}, "
          f"k_e={ke:.4f})")
    print(f"{'k_x':>10}  {'phi_ww(kx,0)':>14}  {'vK shape':>14}  {'ratio':>12}")
    ratios = []
    for kx in kxs:
        num = phi_ww_ky0(kx, p)
        vk = phi_ww_ky0_vK_shape(kx, Lambda)
        ratios.append(num / vk)
        print(f"{kx:10.4f}  {num:14.6e}  {vk:14.6e}  {num / vk:12.6e}")
    ratios = np.array(ratios)
    print(f"ratio constant to rel. spread {np.ptp(ratios) / ratios.mean():.2e}")


if __name__ == "__main__":
    _print_overlay()
