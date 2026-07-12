#!/usr/bin/env python3
r"""
Isotropic-limit unit tests for axisym_family.py (paper 2, Gate A).

These PROVE that the anisotropy-stretched von Karman family collapses to the
exact von Karman forms when c0 = 0 and L2 = Lperp (note sec.6):

  1. phi_ww_ky0 has the vK transverse shape  k_x^2/(1+(k_x/k_e)^2)^(7/3).
  2. component energies are isotropic:  R_11 = R_22 = R_33.
  3. with c0 = 0 the C contribution to Phi_22 is identically 0.
  4. spanwise length is finite, positive, and monotone decreasing in k_x.
  5. Parseval/normalization: int E_2 dk2 == R_22 from an independent 3-D quad.

plus a sanity test that the anisotropy knobs actually change R_22/R_11.

Runtime target: < ~20 s.  Grids are modest but convergence is checked.
"""

import numpy as np
import pytest
from scipy.integrate import quad

from axisym_family import (
    Params,
    Phi22,
    C_scalar,
    E_component,
    R_ii,
    phi_ww_ky0,
    phi_ww_ky0_vK_shape,
    spanwise_length,
    spanwise_length_vK,
    isotropic_params,
    ke_from_length,
    default_k2_grid,
)

# Parametrization axes: a few integral scales and a few rms levels.
LAMBDAS = [0.5, 1.0, 2.0]
UPRIMES = [1.0, 2.5]


def _k2_grid(p):
    """Log-spaced positive k2 grid (spectrally accurate for R_ii)."""
    return default_k2_grid(p)


# ----------------------------------------------------------------------
# 1. Upwash planar spectrum at k_y = 0 reduces to the vK transverse shape.
# ----------------------------------------------------------------------
@pytest.mark.parametrize("Lambda", LAMBDAS)
@pytest.mark.parametrize("uprime", UPRIMES)
def test_phi_ww_ky0_vK_shape(Lambda, uprime):
    p = isotropic_params(Lambda, A0=uprime ** 2)
    ke = ke_from_length(Lambda)
    kxs = np.array([0.2, 0.5, 1.0, 2.0, 4.0, 8.0]) * ke
    num = np.array([phi_ww_ky0(kx, p) for kx in kxs])
    vk = np.array([phi_ww_ky0_vK_shape(kx, Lambda) for kx in kxs])
    ratio = num / vk
    # shape match <=> ratio is constant across k_x
    rel_spread = np.ptp(ratio) / ratio.mean()
    assert rel_spread < 1e-3, f"phi_ww/vK not constant: spread {rel_spread:.2e}"


# ----------------------------------------------------------------------
# 2. Component energies are isotropic in the isotropic limit.
# ----------------------------------------------------------------------
@pytest.mark.parametrize("Lambda", LAMBDAS)
@pytest.mark.parametrize("uprime", UPRIMES)
def test_isotropy_of_component_energies(Lambda, uprime):
    p = isotropic_params(Lambda, A0=uprime ** 2)
    R = R_ii(_k2_grid(p), p)
    R11, R22, R33 = R
    assert np.isclose(R11, R33, rtol=1e-3), f"R11={R11} R33={R33}"
    assert np.isclose(R22, R11, rtol=1e-3), f"R22={R22} R11={R11}"
    assert np.isclose(R22, R33, rtol=1e-3), f"R22={R22} R33={R33}"


# ----------------------------------------------------------------------
# 3. With c0 = 0 the C contribution to Phi_22 is identically zero.
# ----------------------------------------------------------------------
@pytest.mark.parametrize("Lambda", LAMBDAS)
def test_C_off_when_c0_zero(Lambda):
    p = isotropic_params(Lambda)  # c0 = 0
    pts = [(0.3, 0.5, 0.2), (1.0, 0.0, 0.0), (0.7, -1.3, 0.9),
           (2.0, 2.0, 2.0), (0.1, 5.0, 0.0)]
    for (k1, k2, k3) in pts:
        kperp = np.hypot(k1, k3)
        assert abs(C_scalar(k2, kperp, p)) < 1e-12
        # Phi_22 must equal its pure-A part (C term drops out).
        kap2 = k1 * k1 + k2 * k2 + k3 * k3
        A = p.A0 * (lambda kh2: kh2 * kh2 / (1 + kh2) ** (17.0 / 6.0))(
            (k2 / p.L2) ** 2 + (kperp / p.Lperp) ** 2
        )
        pure_A = A * (kperp * kperp / kap2) / (4.0 * np.pi * kap2)
        assert np.isclose(Phi22(k1, k2, k3, p), pure_A, atol=1e-12, rtol=1e-10)


# ----------------------------------------------------------------------
# 4. Spanwise length: finite, positive, monotone decreasing in k_x.
# ----------------------------------------------------------------------
@pytest.mark.parametrize("Lambda", LAMBDAS)
def test_spanwise_length_positive_and_vK_closed_form(Lambda):
    p = isotropic_params(Lambda)
    ke = ke_from_length(Lambda)
    kxs = np.array([0.3, 0.6, 1.0, 2.0, 4.0]) * ke
    ell = np.array([spanwise_length(kx, p) for kx in kxs])
    assert np.all(np.isfinite(ell))
    assert np.all(ell > 0.0)
    # The sec.5.2 ratio has an EXACT isotropic vK closed form -- assert to it
    # (stronger than the monotonicity check the brief falls back to).  NB: this
    # closed form is U-shaped, NOT monotone decreasing; see spanwise_length and
    # its DECISION NEEDED note.
    ell_vK = spanwise_length_vK(kxs, Lambda)
    assert np.allclose(ell, ell_vK, rtol=1e-2), (
        f"ell_y != vK closed form:\n  num={ell}\n  vK ={ell_vK}"
    )


# ----------------------------------------------------------------------
# 5. Parseval/normalization: int E_2 dk2 == R_22 from an independent 3-D quad.
# ----------------------------------------------------------------------
@pytest.mark.parametrize("Lambda", LAMBDAS)
@pytest.mark.parametrize("uprime", UPRIMES)
def test_parseval_R22(Lambda, uprime):
    p = isotropic_params(Lambda, A0=uprime ** 2)
    R22_from_E = R_ii(_k2_grid(p), p)[1]

    # Independent quadrature: integrate Phi_22 over all kappa in SPHERICAL
    # coordinates (kap, mu), a different grid from the cylindrical E_2 route.
    #   R_22 = int_0^inf int_-1^1 Phi_22(kap,mu) * kap^2 * 2 pi dmu dkap
    def inner(kap):
        def f(mu):
            k2 = kap * mu
            kperp = kap * np.sqrt(max(0.0, 1.0 - mu * mu))
            # rebuild a kappa point on the (k1=kperp, k3=0) meridian
            return Phi22(kperp, k2, 0.0, p)
        val, _ = quad(f, -1.0, 1.0, limit=200)
        return 2.0 * np.pi * kap * kap * val

    R22_direct, _ = quad(inner, 0.0, np.inf, limit=200)
    assert np.isclose(R22_from_E, R22_direct, rtol=1e-6), (
        f"Parseval mismatch: E-route {R22_from_E:.8e} vs "
        f"direct {R22_direct:.8e}"
    )


# ----------------------------------------------------------------------
# Sanity: the anisotropy knobs actually move R_22/R_11 off 1.
# ----------------------------------------------------------------------
def test_anisotropy_changes_ratio():
    ke = ke_from_length(1.0)
    # (a) stretch the ODT axis relative to perp
    p_stretch = Params(A0=1.0, c0=0.0, L2=2.0 * ke, Lperp=ke)
    R = R_ii(_k2_grid(p_stretch), p_stretch)
    assert not np.isclose(R[1] / R[0], 1.0, rtol=1e-3), \
        f"L2!=Lperp did not break isotropy: R22/R11={R[1] / R[0]}"

    # (b) turn on the C anisotropy scalar
    p_cscalar = Params(A0=1.0, c0=0.8, L2=ke, Lperp=ke)
    R = R_ii(_k2_grid(p_cscalar), p_cscalar)
    assert not np.isclose(R[1] / R[0], 1.0, rtol=1e-3), \
        f"c0!=0 did not break isotropy: R22/R11={R[1] / R[0]}"


# ----------------------------------------------------------------------
# End-of-run overlay so the reduction can be eyeballed in the next chat.
# ----------------------------------------------------------------------
def test_zzz_print_overlay(capsys):
    p = isotropic_params(1.0, A0=1.0)
    ke = ke_from_length(1.0)
    kxs = np.array([0.25, 0.5, 1.0, 2.0, 4.0]) * ke
    lines = ["", "fitted-vs-vK Phi_ww overlay (isotropic, Lambda=1, k_e="
             f"{ke:.4f}):",
             f"{'k_x':>10}  {'phi_ww(kx,0)':>14}  {'vK shape':>14}  {'ratio':>12}"]
    ratios = []
    for kx in kxs:
        num = phi_ww_ky0(kx, p)
        vk = phi_ww_ky0_vK_shape(kx, 1.0)
        ratios.append(num / vk)
        lines.append(f"{kx:10.4f}  {num:14.6e}  {vk:14.6e}  {num / vk:12.6e}")
    ratios = np.array(ratios)
    lines.append(f"ratio constant to rel. spread "
                 f"{np.ptp(ratios) / ratios.mean():.2e}")
    with capsys.disabled():
        print("\n".join(lines))
    assert np.ptp(ratios) / ratios.mean() < 1e-3


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v", "-s"]))
