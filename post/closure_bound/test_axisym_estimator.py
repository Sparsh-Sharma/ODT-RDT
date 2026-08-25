"""Isotropic-limit test suite for the LOS closure-bound estimator.

Every assertion here is an exact statement from the manuscript (Section 5) or
from classical isotropic-turbulence kinematics; tolerances reflect quadrature
error only.
"""

import numpy as np
import pytest

from axisym_estimator import (
    Grid,
    a_isotropic,
    bound_pi,
    g_kernels_a,
    g_kernels_c,
    isotropic_calibration,
    kinetic_energy,
    line_spectra,
    pi_exact,
)


# ----------------------------------------------------------------------------
# kernel identities (exact, manuscript eqs. 5.8-5.10 and 5.11-5.12)
# ----------------------------------------------------------------------------

X = np.linspace(0.0, 1.0, 501)


def test_kernels_traceless_pointwise():
    for gk in (g_kernels_a, g_kernels_c):
        g1, g2, g3 = gk(X)
        np.testing.assert_allclose(g1 + g2 + g3, 0.0, atol=1e-14)


def test_kernels_vanish_at_aligned_wavevector():
    for gk in (g_kernels_a, g_kernels_c):
        for g in gk(np.array([1.0])):
            assert abs(g[0]) < 1e-14


def test_kernel_isotropic_angular_averages():
    # <f(x)> over the sphere = int_0^1 f(mu^2) dmu ; expect +1/5 : -1/5 : 0.
    mu = np.linspace(0.0, 1.0, 200001)
    g1, g2, g3 = g_kernels_a(mu**2)
    avg = [np.trapezoid(g, mu) for g in (g1, g2, g3)]
    np.testing.assert_allclose(avg, [0.2, -0.2, 0.0], atol=1e-9)


# ----------------------------------------------------------------------------
# forward map, isotropic limit
# ----------------------------------------------------------------------------

@pytest.fixture(scope="module")
def iso():
    grid = Grid.make(k_min=1e-3, k_max=1e3, n_k2=160, n_kperp=320)
    a_fn = lambda k2v, kp: a_isotropic(k2v, kp, L=1.0, amplitude=1.0)
    phi11, phi22 = line_spectra(grid, a_fn)
    return grid, a_fn, phi11, phi22


def test_forward_map_energy_normalization(iso):
    # int (phi_11 + phi_22 + phi_33) dk2 over full k2 axis = R_mm = 2 k_t.
    grid, a_fn, phi11, phi22 = iso
    r_mm = 2.0 * np.sum(grid.w2 * (2.0 * phi11 + phi22))  # 2x: k2<0 half
    kt = kinetic_energy(grid, a_fn)
    np.testing.assert_allclose(r_mm, 2.0 * kt, rtol=1e-6)


def test_forward_map_isotropic_derivative_relation(iso):
    # Classical isotropic relation between transverse and longitudinal 1-D
    # spectra along the line: phi_T = (1/2) (phi_L - k2 dphi_L/dk2),
    # with phi_L = phi_22 (u_2 is longitudinal on the x_2 line).
    grid, _, phi11, phi22 = iso
    # differentiate in ln(k2): uniform spacing on the geomspace grid keeps
    # the central difference genuinely second-order accurate
    dphi22 = np.gradient(phi22, np.log(grid.k2)) / grid.k2
    predicted_t = 0.5 * (phi22 - grid.k2 * dphi22)
    sl = slice(20, -20)  # avoid grid-edge quadrature/gradient artifacts
    # residual ~2e-3 is trapezoid quadrature in the vK tail (flat in k2,
    # insensitive to differencing); a wrong angular weight would be O(1)
    np.testing.assert_allclose(phi11[sl], predicted_t[sl], rtol=5e-3)


def test_forward_map_isotropy_of_transverse(iso):
    # phi_11 = phi_33 holds by construction under A1; here confirm both
    # exceed 0 and that longitudinal < transverse at high k2 (classical).
    grid, _, phi11, phi22 = iso
    assert np.all(phi11 > 0.0) and np.all(phi22 > 0.0)
    high = grid.k2 > 10.0
    assert np.all(phi11[high] > phi22[high])


# ----------------------------------------------------------------------------
# exact rapid pressure-strain, isotropic limit (Crow 1968)
# ----------------------------------------------------------------------------

def test_exact_pi_isotropic(iso):
    # Pi = (4/5) k_t S_ij with S = diag(1/2,-1/2,0) -> (2/5, -2/5, 0) k_t.
    grid, a_fn, _, _ = iso
    kt = kinetic_energy(grid, a_fn)
    pi = pi_exact(grid, a_fn)
    np.testing.assert_allclose(pi / kt, [0.4, -0.4, 0.0], atol=2e-3)


def test_exact_pi_traceless_with_anisotropy(iso):
    # tracelessness must survive an arbitrary realizable c-scalar.
    grid, a_fn, _, _ = iso
    c_fn = lambda k2v, kp: 0.5 * a_isotropic(k2v, kp)  # a + (1-x)c > 0 ok
    pi = pi_exact(grid, a_fn, c_fn)
    kt = kinetic_energy(grid, a_fn, c_fn)
    assert abs(pi.sum()) / kt < 1e-10


# ----------------------------------------------------------------------------
# the LP bounds (the estimator itself)
# ----------------------------------------------------------------------------

@pytest.fixture(scope="module")
def lp_iso():
    # coarser grid: one LP pair per (k2, component); keep runtime modest
    grid = Grid.make(k_min=1e-2, k_max=1e2, n_k2=40, n_kperp=120)
    a_fn = lambda k2v, kp: a_isotropic(k2v, kp, L=1.0, amplitude=1.0)
    phi11, phi22 = line_spectra(grid, a_fn)
    bounds = bound_pi(grid, phi11, phi22)
    return grid, a_fn, bounds


def test_lp_all_solved(lp_iso):
    _, _, bounds = lp_iso
    assert bounds.status_ok.all()


def test_lp_bounds_ordered(lp_iso):
    _, _, bounds = lp_iso
    assert np.all(bounds.upper >= bounds.lower - 1e-12)


def test_lp_bounds_bracket_truth_pointwise(lp_iso):
    # the true isotropic spectrum is feasible, so at every k2 the exact
    # pi_n density must lie inside [lower, upper].
    grid, a_fn, bounds = lp_iso
    from axisym_estimator import PI_PREFACTOR, g_kernels_a as gka, \
        g_kernels_c as gkc  # noqa: N813

    for i, k2v in enumerate(grid.k2):
        x = grid.x(k2v)
        a = a_fn(k2v, grid.kperp)
        for comp in range(3):
            density = PI_PREFACTOR * np.sum(grid.wperp * gka(x)[comp] * a)
            assert bounds.lower[comp, i] - 1e-9 <= density
            assert density <= bounds.upper[comp, i] + 1e-9


def test_lp_integrated_bounds_bracket_crow(lp_iso):
    grid, a_fn, bounds = lp_iso
    kt = kinetic_energy(grid, a_fn)
    lo, hi = bounds.integrated(grid.w2)
    target = np.array([0.4, -0.4, 0.0]) * kt
    tol = 5e-3 * kt
    assert np.all(lo <= target + tol) and np.all(target - tol <= hi)


def test_lp_bounds_scale_linearly(lp_iso):
    # doubling the spectrum amplitude must double the bounds (linearity).
    grid, _, bounds = lp_iso
    a2 = lambda k2v, kp: a_isotropic(k2v, kp, L=1.0, amplitude=2.0)
    phi11, phi22 = line_spectra(grid, a2)
    b2 = bound_pi(grid, phi11, phi22)
    np.testing.assert_allclose(b2.lower, 2.0 * bounds.lower,
                               rtol=1e-6, atol=1e-12)
    np.testing.assert_allclose(b2.upper, 2.0 * bounds.upper,
                               rtol=1e-6, atol=1e-12)


def test_lp_band_widens_bounds(lp_iso):
    # replacing equalities by a data band can only widen the interval.
    grid, a_fn, bounds = lp_iso
    phi11, phi22 = line_spectra(grid, a_fn)
    banded = bound_pi(grid, phi11, phi22, data_band=0.05)
    assert np.all(banded.lower <= bounds.lower + 1e-9)
    assert np.all(banded.upper >= bounds.upper - 1e-9)


def test_anisotropy_moves_line_data(iso):
    # sanity: a realizable c-scalar changes phi_22/phi_11 (the knob works).
    grid, a_fn, phi11, phi22 = iso
    c_fn = lambda k2v, kp: 0.5 * a_isotropic(k2v, kp)
    p11c, p22c = line_spectra(grid, a_fn, c_fn)
    ratio_iso = phi22 / phi11
    ratio_ani = p22c / p11c
    assert np.max(np.abs(ratio_ani - ratio_iso)) > 0.05


# ----------------------------------------------------------------------------
# stage b: log-Lipschitz slope cap
# ----------------------------------------------------------------------------

@pytest.fixture(scope="module")
def lp_grid():
    grid = Grid.make(k_min=1e-2, k_max=1e2, n_k2=40, n_kperp=120)
    a_fn = lambda k2v, kp: a_isotropic(k2v, kp, L=1.0, amplitude=1.0)
    phi11, phi22 = line_spectra(grid, a_fn)
    return grid, a_fn, phi11, phi22


def test_stage_b_none_reproduces_stage_a(lp_iso, lp_grid):
    # slope_cap=None must give the raw (stage-a) bounds bit-for-bit in value
    # despite the internal c = cp - cm split.
    _, _, raw = lp_iso
    grid, _, phi11, phi22 = lp_grid
    again = bound_pi(grid, phi11, phi22, slope_cap=None)
    np.testing.assert_allclose(again.lower, raw.lower, rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(again.upper, raw.upper, rtol=1e-8, atol=1e-12)


def test_stage_b_truth_still_feasible(lp_grid):
    # with lambda = 4 > 11/3 the vK truth satisfies the cap, so the bounds
    # must still bracket the exact pi_n density at every k2.
    grid, a_fn, phi11, phi22 = lp_grid
    from axisym_estimator import PI_PREFACTOR, g_kernels_a as gka

    bounds = bound_pi(grid, phi11, phi22, slope_cap=4.0)
    assert bounds.status_ok.all()
    for i, k2v in enumerate(grid.k2):
        x = grid.x(k2v)
        a = a_fn(k2v, grid.kperp)
        for comp in range(3):
            density = PI_PREFACTOR * np.sum(grid.wperp * gka(x)[comp] * a)
            assert bounds.lower[comp, i] - 1e-9 <= density
            assert density <= bounds.upper[comp, i] + 1e-9


def test_stage_b_tightens_monotonically(lp_grid):
    # smaller lambda = stronger constraint = narrower (nested) bounds,
    # up to LP solve tolerance (relative).
    grid, _, phi11, phi22 = lp_grid
    b_raw = bound_pi(grid, phi11, phi22)
    b8 = bound_pi(grid, phi11, phi22, slope_cap=8.0)
    b4 = bound_pi(grid, phi11, phi22, slope_cap=4.0)
    for tight, loose in ((b8, b_raw), (b4, b8)):
        tol = 1e-6 * np.abs(loose.upper - loose.lower).max()
        assert np.all(tight.lower >= loose.lower - tol)
        assert np.all(tight.upper <= loose.upper + tol)


def test_stage_b_slope_cap_alone_is_weakly_binding(lp_grid):
    # FINDING (2026-08-25, isotropic calibration): a log-Lipschitz cap in
    # kperp, even at the tightest vK-feasible lambda=4, narrows the
    # integrated width by <5%.  A log-Lipschitz spectrum can still vary
    # like kperp^{+-4}, so angular reweighting survives smoothing.  This
    # test pins the finding; if it ever fails, the conclusion in the note
    # must be revisited.
    grid, _, phi11, phi22 = lp_grid
    lo_r, hi_r = bound_pi(grid, phi11, phi22).integrated(grid.w2)
    lo_4, hi_4 = bound_pi(grid, phi11, phi22,
                          slope_cap=4.0).integrated(grid.w2)
    assert np.all((hi_4 - lo_4) > 0.95 * (hi_r - lo_r))


def test_stage_b_polarization_cap_is_strongly_binding(lp_grid):
    # FINDING: the polarization cap |c| <= gamma a is the lever.  gamma=0
    # (isotropic polarization; exact for the isotropic truth, so bounds
    # must still bracket) cuts every integrated width by >40%.
    grid, a_fn, phi11, phi22 = lp_grid
    from axisym_estimator import PI_PREFACTOR, g_kernels_a as gka

    b_raw = bound_pi(grid, phi11, phi22)
    b_pol = bound_pi(grid, phi11, phi22, polarization_cap=0.0)
    assert b_pol.status_ok.all()
    # truth (c=0) satisfies the cap -> still bracketed pointwise
    for i, k2v in enumerate(grid.k2):
        x = grid.x(k2v)
        a = a_fn(k2v, grid.kperp)
        for comp in range(3):
            density = PI_PREFACTOR * np.sum(grid.wperp * gka(x)[comp] * a)
            assert b_pol.lower[comp, i] - 1e-9 <= density
            assert density <= b_pol.upper[comp, i] + 1e-9
    lo_r, hi_r = b_raw.integrated(grid.w2)
    lo_p, hi_p = b_pol.integrated(grid.w2)
    assert np.all((hi_p - lo_p) < 0.6 * (hi_r - lo_r))


# ----------------------------------------------------------------------------
# the calibration number
# ----------------------------------------------------------------------------

def test_isotropic_calibration_report():
    cal = isotropic_calibration(n_k2=40, n_kperp=120,
                                k_min=1e-2, k_max=1e2)
    kt = cal["k_t"]
    print("\nIsotropic vK calibration (the HIT-limit indeterminacy):")
    print("comp   exact/kt   lower/kt   upper/kt     I_n")
    for n, name in enumerate(("11", "22", "33")):
        print(f"Pi_{name}  {cal['pi_exact'][n]/kt:+8.4f}  "
              f"{cal['pi_lower'][n]/kt:+8.4f}  "
              f"{cal['pi_upper'][n]/kt:+8.4f}  {cal['indeterminacy'][n]:8.4f}")
    assert np.all(cal["indeterminacy"] >= 0.0)
    assert np.all(np.isfinite(cal["indeterminacy"]))
