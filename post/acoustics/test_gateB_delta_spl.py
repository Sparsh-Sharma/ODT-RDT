#!/usr/bin/env python3
r"""
Tests for gateB_delta_spl.py: Amiet spanwise length closed form, Delta-SPL
identities (self-baseline zero, amplitude level shift), and the anisotropy
spectral-shape signature (paper 2, Session 3).

No run data is read.  Fast.
"""

import numpy as np
import pytest

import axisym_family as af
from axisym_family import Params
import gateB_delta_spl as g

LAMBDAS = [0.5, 1.0, 2.0]
U = 30.0


# ----------------------------------------------------------------------
# 1. Amiet spanwise correlation length matches the exact isotropic vK form.
# ----------------------------------------------------------------------
@pytest.mark.parametrize("Lambda", LAMBDAS)
def test_spanwise_corr_length_vK(Lambda):
    p = af.isotropic_params(Lambda, A0=1.3)
    ke = af.ke_from_length(Lambda)
    for r in [0.3, 0.6, 1.0, 2.0, 4.0, 8.0]:
        kx = r * ke
        num = g.spanwise_corr_length(kx, p)
        cf = g.spanwise_corr_length_vK(kx, Lambda)
        assert num > 0.0
        assert np.isclose(num, cf, rtol=1e-3), (r, num, cf)


# ----------------------------------------------------------------------
# 2. High-frequency regime: ell_y decreases (like 1/k_x) for k_x >= k_e.
# ----------------------------------------------------------------------
@pytest.mark.parametrize("Lambda", LAMBDAS)
def test_spanwise_length_high_freq_decay(Lambda):
    p = af.isotropic_params(Lambda)
    ke = af.ke_from_length(Lambda)
    kxs = np.array([1.0, 1.5, 2.0, 3.0, 5.0, 8.0]) * ke
    ell = np.array([g.spanwise_corr_length(kx, p) for kx in kxs])
    assert np.all(np.diff(ell) < 0.0), ell
    # asymptotic ~1/k_x: ell_y(2 k_x)/ell_y(k_x) -> 1/2 at large k_x
    hi = 20.0 * ke
    ratio = g.spanwise_corr_length(2 * hi, p) / g.spanwise_corr_length(hi, p)
    assert np.isclose(ratio, 0.5, rtol=5e-2), ratio


# ----------------------------------------------------------------------
# 3. Delta-SPL against the same field is identically zero.
# ----------------------------------------------------------------------
def test_delta_spl_self_baseline_zero():
    p = af.isotropic_params(1.0, A0=1.7)
    omega = np.geomspace(1.0, 500.0, 12)
    d = g.delta_spl(omega, U, p, p)
    assert np.allclose(d, 0.0, atol=1e-9)


# ----------------------------------------------------------------------
# 4. Pure amplitude change is a flat level shift 10 log10(A0_a/A0_i).
# ----------------------------------------------------------------------
@pytest.mark.parametrize("alpha", [2.0, 0.5, 3.3])
def test_delta_spl_amplitude_level_shift(alpha):
    p_i = af.isotropic_params(1.0, A0=1.0)
    p_a = af.isotropic_params(1.0, A0=alpha)      # only A0 differs
    omega = np.geomspace(1.0, 500.0, 10)
    d = g.delta_spl(omega, U, p_a, p_i)
    assert np.allclose(d, 10.0 * np.log10(alpha), atol=1e-6)


# ----------------------------------------------------------------------
# 5. noise kernel is linear in A0 (energy).
# ----------------------------------------------------------------------
def test_noise_kernel_linear_in_A0():
    ke = af.ke_from_length(1.0)
    p1 = Params(A0=1.0, c0=0.0, L2=ke, Lperp=ke)
    p2 = Params(A0=2.5, c0=0.0, L2=ke, Lperp=ke)
    for kx in [0.3 * ke, ke, 3 * ke]:
        assert np.isclose(g.noise_kernel(kx, p2) / g.noise_kernel(kx, p1), 2.5,
                          rtol=1e-8)


# ----------------------------------------------------------------------
# 6. Anisotropy is a frequency-dependent (spectral-shape) signature, unlike a
#    pure amplitude shift.
# ----------------------------------------------------------------------
def test_anisotropy_is_frequency_dependent():
    ke = af.ke_from_length(1.0)
    p_a = Params(A0=1.0, c0=0.0, L2=0.5 * ke, Lperp=ke)   # stretched k2 scale
    p_i = g.matched_isotropic(p_a, mode="perp")
    omega = np.geomspace(0.2, 40.0, 9) * U * ke
    d = g.delta_spl(omega, U, p_a, p_i)
    assert np.all(np.isfinite(d))
    # not a constant level shift: low- and high-frequency ends differ clearly
    assert abs(d[0] - d[-1]) > 1.0, d
    # a c0-driven anisotropy is likewise frequency dependent and finite
    p_c = Params(A0=1.0, c0=0.8, L2=ke, Lperp=ke)
    d_c = g.delta_spl(omega, U, p_c, g.matched_isotropic(p_c, mode="perp"))
    assert np.all(np.isfinite(d_c)) and np.ptp(d_c) > 0.0


# ----------------------------------------------------------------------
# 7. matched_isotropic really is isotropic.
# ----------------------------------------------------------------------
@pytest.mark.parametrize("mode", ["perp", "mean"])
def test_matched_isotropic_is_isotropic(mode):
    p = Params(A0=1.4, c0=0.7, L2=0.6, Lperp=1.5)
    b = g.matched_isotropic(p, mode=mode)
    assert b.c0 == 0.0
    assert np.isclose(b.L2, b.Lperp)
    assert b.A0 == p.A0


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
