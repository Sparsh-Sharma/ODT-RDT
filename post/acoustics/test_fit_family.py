#!/usr/bin/env python3
r"""
Synthetic self-consistency tests for fit_family.py (paper 2, Gate A, Session 2).

These PROVE the fitter can recover the four family parameters from the line
spectra they generate -- the identifiability/round-trip check that de-risks the
fit before any real ODT run data is attached:

  1. noise-free anisotropic recovery: fit(E_i(p*)) -> p*  to tight tolerance.
  2. isotropic data recovers isotropy: L2 ~ Lperp and c0 ~ 0.
  3. robustness to moderate spectral noise (recovery within a looser band).
  4. forward-map sanity: E1 == E3, and (1/2) sum R_ii = kt.

No run data is read.  Runtime target: a few seconds.
"""

import numpy as np
import pytest

import axisym_family as af
from axisym_family import Params
from fit_family import (
    line_spectra,
    observables,
    fit_family,
    synthetic_targets,
)


def _grid(p, half=18, n=140, thin=4):
    return af.default_k2_grid(p, half=half, n=n)[::thin]


def _rel(a, b):
    denom = abs(a) if abs(a) > 1e-9 else 1.0
    return abs(b - a) / denom


# ----------------------------------------------------------------------
# 1. Noise-free recovery of known anisotropic parameters.
# ----------------------------------------------------------------------
@pytest.mark.parametrize("p_true", [
    Params(A0=1.0, c0=0.0, L2=1.0, Lperp=1.0),
    Params(A0=2.3, c0=0.0, L2=0.7, Lperp=1.6),   # length anisotropy only
    Params(A0=1.4, c0=0.6, L2=0.8, Lperp=1.3),   # length + c0 anisotropy
])
def test_recover_params_noise_free(p_true):
    k2 = _grid(p_true)
    E1, E2 = synthetic_targets(p_true, k2, noise=0.0)
    res = fit_family(k2, E1, E2)
    assert res.success
    assert _rel(p_true.A0, res.params.A0) < 1e-2
    assert _rel(p_true.L2, res.params.L2) < 1e-2
    assert _rel(p_true.Lperp, res.params.Lperp) < 1e-2
    # c0 enters only through the (small) anisotropy correction -> looser abs tol
    assert abs(res.params.c0 - p_true.c0) < 2e-2


# ----------------------------------------------------------------------
# 2. Isotropic data must recover isotropy (L2 ~ Lperp, c0 ~ 0) from an
#    anisotropic starting guess.
# ----------------------------------------------------------------------
def test_isotropic_data_recovers_isotropy():
    Lam = 1.0
    p_true = af.isotropic_params(Lam, A0=1.7)
    k2 = _grid(p_true)
    E1, E2 = synthetic_targets(p_true, k2, noise=0.0)
    # deliberately anisotropic, off-amplitude start
    p0 = Params(A0=0.5, c0=0.5, L2=0.5 * p_true.L2, Lperp=2.0 * p_true.Lperp)
    res = fit_family(k2, E1, E2, p0=p0)
    assert res.success
    assert _rel(res.params.L2, res.params.Lperp) < 1e-2
    assert abs(res.params.c0) < 2e-2


# ----------------------------------------------------------------------
# 3. Robustness to moderate multiplicative spectral noise.
# ----------------------------------------------------------------------
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_recovery_under_noise(seed):
    p_true = Params(A0=1.5, c0=0.0, L2=0.9, Lperp=1.4)
    k2 = _grid(p_true, thin=3)
    E1, E2 = synthetic_targets(p_true, k2, noise=0.05, seed=seed)
    res = fit_family(k2, E1, E2)
    assert res.success
    # 5% spectral noise -> length scales still within a few percent
    assert _rel(p_true.L2, res.params.L2) < 5e-2
    assert _rel(p_true.Lperp, res.params.Lperp) < 5e-2
    assert _rel(p_true.A0, res.params.A0) < 8e-2


# ----------------------------------------------------------------------
# 4. Forward-map sanity.
# ----------------------------------------------------------------------
def test_forward_map_consistency():
    p = Params(A0=1.3, c0=0.4, L2=0.9, Lperp=1.2)
    k2 = af.default_k2_grid(p)
    obs = observables(p, k2)
    assert np.allclose(obs["E1"], obs["E3"])          # axisymmetry
    assert np.isclose(obs["kt"], 0.5 * np.sum(obs["R"]))
    E1, E2 = line_spectra(p, k2[::20])
    assert np.all(E1 > 0) and np.all(E2 > 0)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
