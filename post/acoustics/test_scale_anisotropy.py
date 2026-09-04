#!/usr/bin/env python3
r"""
Tests for scale_anisotropy.py (Test-3 rate-competition diagnostic).

All synthetic; no run data.  Covers:
  1. the isotropic reference ratio's exact limits (2 at low k, 3/4 at high k),
  2. isotropic data -> A ~= 1 everywhere and a gracefully degenerate decay fit,
  3. recovery of a controlled, imposed anisotropy-decay exponent,
  4. scale-persistent (family) anisotropy -> slope ~= 0 ("Test-3 signature"),
  5. end-to-end: dump file -> scale_anisotropy() recovers the imposed decay.
"""

import os

import numpy as np
import pytest

import axisym_family as af
from axisym_family import Params
import odt_io
import scale_anisotropy as sa

KE = af.ke_from_length(1.0)


def _iso_spectra(k2, half=22, n=400):
    p = Params(A0=1.0, c0=0.0, L2=KE, Lperp=KE)
    E1 = af.E_component(1, k2, p, half=half, n=n)
    E2 = af.E_component(2, k2, p, half=half, n=n)
    return E1, E2


# ----------------------------------------------------------------------
# 1. Isotropic reference ratio limits: 2 (low k) and 3/4 (the 4/3 law).
# ----------------------------------------------------------------------
def test_iso_reference_limits():
    lo = sa.iso_reference_ratio(np.array([1e-3 * KE]), KE)
    hi = sa.iso_reference_ratio(np.array([1e3 * KE]), KE)
    assert np.isclose(lo.item(), 2.0, rtol=1e-3)
    assert np.isclose(hi.item(), 0.75, rtol=1e-3)


# ----------------------------------------------------------------------
# 2. Isotropic data: A ~= 1, decay fit degenerates gracefully.
# ----------------------------------------------------------------------
def test_isotropic_data_gives_unity():
    k2 = af.default_k2_grid(Params(1, 0, KE, KE), half=14, n=100)[::2]
    E1, E2 = _iso_spectra(k2)
    an = sa.anisotropy_function(k2, E1, E2, E1)
    assert np.allclose(an["A"], 1.0, atol=2e-3)
    assert np.isclose(an["L_iso"], KE, rtol=1e-2)
    dc = sa.decay_per_octave(k2, an["A"], floor=1e-2)
    assert dc.npts < 4 and np.isnan(dc.slope)   # nothing above the floor


# ----------------------------------------------------------------------
# 3. Controlled decay: impose |A-1| ~ k^-r on top of isotropy, recover r.
# ----------------------------------------------------------------------
@pytest.mark.parametrize("r", [0.3, 0.5, 1.0])
def test_recover_imposed_decay(r):
    k2 = af.default_k2_grid(Params(1, 0, KE, KE), half=14, n=160)[::2]
    E1, E2i = _iso_spectra(k2)
    E2 = E2i * (1.0 + 0.4 * (k2 / KE) ** (-r))
    an = sa.anisotropy_function(k2, E1, E2, E1)     # L_iso fitted, not given
    dc = sa.decay_per_octave(k2, an["A"], kband=(2 * KE, 50 * KE))
    assert dc.npts >= 6
    assert np.isclose(dc.slope, -r, atol=0.03), dc
    assert np.isclose(dc.factor_per_map, 3.0 ** (-r), rtol=0.05)


# ----------------------------------------------------------------------
# 4. Scale-persistent anisotropy (the family itself): slope ~= 0.
# ----------------------------------------------------------------------
def test_persistent_anisotropy_slope_zero():
    p = Params(A0=1.0, c0=0.0, L2=0.6 * KE, Lperp=1.3 * KE)  # pure stretch
    k2 = af.default_k2_grid(p, half=14, n=160)[::2]
    E1 = af.E_component(1, k2, p, half=22, n=400)
    E2 = af.E_component(2, k2, p, half=22, n=400)
    an = sa.anisotropy_function(k2, E1, E2, E1)
    dc = sa.decay_per_octave(k2, an["A"], kband=(3 * KE, 60 * KE))
    assert dc.npts >= 6
    assert abs(dc.slope) < 0.08, dc                  # persists across scales
    band = (k2 >= 3 * KE) & (k2 <= 60 * KE)
    assert np.all(np.abs(an["A"][band] - 1.0) > 0.2)  # and is clearly nonzero


# ----------------------------------------------------------------------
# 5. End to end: dump -> scale_anisotropy recovers the imposed decay.
# ----------------------------------------------------------------------
def _write_dump_with_spectra(fn, Ldom, N, E1f, E2f, E3f, seed=0):
    """Uniform-grid dump whose component spectral densities are E{1,2,3}f(k)."""
    rng = np.random.default_rng(seed)
    faces = np.linspace(-0.5 * Ldom, 0.5 * Ldom, N + 1)
    posf = faces[:-1]
    kfft = 2.0 * np.pi * np.fft.rfftfreq(N, d=Ldom / N)
    dk = 2.0 * np.pi / Ldom
    fields = []
    for Ef in (E1f, E2f, E3f):
        amp = np.zeros_like(kfft)
        amp[1:] = np.sqrt(np.maximum(Ef(kfft[1:]), 0.0) * dk / 2.0)
        phase = np.exp(2j * np.pi * rng.random(kfft.size))
        fields.append(np.fft.irfft(amp * phase, n=N) * N)
    odt_io.write_dump(fn, posf, *fields, t=1.0)


def test_end_to_end_dump_driver(tmp_path):
    r, d0 = 0.5, 0.4
    Ldom, N = 40.0, 4096
    p_iso = Params(A0=1.0, c0=0.0, L2=KE, Lperp=KE)
    E1f = lambda k: af.E_component(1, k, p_iso, half=22, n=400)
    E2f = lambda k: (af.E_component(2, k, p_iso, half=22, n=400)
                     * (1.0 + d0 * (k / KE) ** (-r)))
    fn = os.path.join(tmp_path, "dmp_0001.dat")
    _write_dump_with_spectra(fn, Ldom, N, E1f, E2f, E1f)

    res = sa.scale_anisotropy(str(tmp_path), Nu=N, nbin=60,
                              kmin=0.3 * KE, kmax=60 * KE,
                              kband=(2 * KE, 40 * KE))
    dc = res["decay"]
    assert dc.npts >= 8
    assert np.isclose(dc.slope, -r, atol=0.06), dc
    assert np.isclose(dc.factor_per_map, 3.0 ** (-r), rtol=0.10)
    assert isinstance(sa.report(res), str) and "per triplet map" in sa.report(res)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))


# ----------------------------------------------------------------------
# 6. ref="equal": component-equal spectra are exactly isotropic (A = 1).
# ----------------------------------------------------------------------
def test_equal_reference():
    k2 = af.default_k2_grid(Params(1, 0, KE, KE), half=12, n=60)
    E1, E2 = _iso_spectra(k2)
    an = sa.anisotropy_function(k2, E1, E1, E1, ref="equal")
    assert np.allclose(an["A"], 1.0)
    an2 = sa.anisotropy_function(k2, E1, 2.0 * E1, E1, ref="equal")
    assert np.allclose(an2["A"], 2.0)      # rho_iso = 1: A is the raw ratio
    with pytest.raises(ValueError):
        sa.anisotropy_function(k2, E1, E2, E1, ref="bogus")
