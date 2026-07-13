#!/usr/bin/env python3
r"""
Tests for odt_io.py: dump round-trip, spectral-density normalization (Parseval),
single-mode recovery, and the fit-target loader wiring (paper 2, Session 2).

No real run data is read; synthetic dumps are written by odt_io.write_dump.
"""

import os

import numpy as np
import pytest

import odt_io


def _uniform_line(N=1024, L=2.0):
    """Face positions for a domain centred at 0 with length L (posf[0] = -L/2)."""
    faces = np.linspace(-0.5 * L, 0.5 * L, N + 1)
    posf = faces[:-1]                       # read_dump/spectra use posf[0]=-L/2
    xc = 0.5 * (faces[:-1] + faces[1:])
    return posf, xc, L


def test_dump_round_trip(tmp_path):
    posf, xc, L = _uniform_line(64)
    u = np.sin(xc)
    v = np.cos(2 * xc)
    w = 0.3 * xc
    fn = os.path.join(tmp_path, "dmp_0001.dat")
    odt_io.write_dump(fn, posf, u, v, w, t=1.234)
    t, p2, u2, v2, w2 = odt_io.read_dump(fn)
    assert np.isclose(t, 1.234)
    assert np.allclose(p2, posf)
    assert np.allclose(u2, u) and np.allclose(v2, v) and np.allclose(w2, w)


def test_single_mode_recovery(tmp_path):
    # a pure Fourier mode of wavenumber k0 along the line -> spectral peak at k0
    N, L = 2048, 2.0 * np.pi        # so fundamental dk = 2*pi/L = 1
    posf, xc, L = _uniform_line(N, L)
    m = 5
    k0 = 2.0 * np.pi * m / L
    u = np.zeros_like(xc)
    v = np.sqrt(2.0) * np.cos(k0 * xc)     # amplitude -> variance 1
    w = np.zeros_like(xc)
    fn = os.path.join(tmp_path, "dmp_0002.dat")
    odt_io.write_dump(fn, posf, u, v, w, t=0.0)
    _, posf2, u2, v2, w2 = odt_io.read_dump(fn)
    k2, E1, E2, E3 = odt_io.component_spectra(posf2, u2, v2, w2, Nu=N)
    assert np.isclose(k2[np.argmax(E2)], k0, rtol=2e-2)
    # u,w channels are silent
    assert E1.max() < 1e-6 * E2.max()
    assert E3.max() < 1e-6 * E2.max()


def test_density_parseval(tmp_path):
    # broadband random field: int E_i dk2 ~= <u_i^2>
    N, L = 4096, 3.0
    posf, xc, L = _uniform_line(N, L)
    rng = np.random.default_rng(0)
    u = rng.normal(0, 1.0, N)
    v = rng.normal(0, 0.5, N)
    w = rng.normal(0, 2.0, N)
    fn = os.path.join(tmp_path, "dmp_0003.dat")
    odt_io.write_dump(fn, posf, u, v, w, t=0.0)
    _, p2, u2, v2, w2 = odt_io.read_dump(fn)
    k2, E1, E2, E3 = odt_io.component_spectra(p2, u2, v2, w2, Nu=N)
    # trapz of density over k adds back the k=0 bin's absence negligibly
    for E, f in ((E1, u2), (E2, v2), (E3, w2)):
        integ = np.trapezoid(E, k2)
        var = np.var(f - f.mean())
        assert np.isclose(integ, var, rtol=0.15), (integ, var)


def test_load_fit_target_and_band(tmp_path):
    N, L = 2048, 2.0
    posf, xc, L = _uniform_line(N, L)
    rng = np.random.default_rng(1)
    u = rng.normal(0, 1, N); v = rng.normal(0, 1, N); w = rng.normal(0, 1, N)
    for i, t in enumerate([0.0, 1.0, 2.0]):
        odt_io.write_dump(os.path.join(tmp_path, f"dmp_{i:04d}.dat"),
                          posf, u, v, w, t=t)
    d = odt_io.load_fit_target(str(tmp_path), strain=1.0, smag=1.0, Nu=N)
    assert np.isclose(d["e"], 1.0) and np.isclose(d["L"], L)
    assert d["k2"].size == d["E2"].size > 0

    k2, Eperp, E2 = odt_io.target_for_fit(str(tmp_path), strain=1.0, Nu=N,
                                          kmin=2.0, kmax=50.0)
    assert k2.min() >= 2.0 and k2.max() <= 50.0
    assert Eperp.shape == E2.shape == k2.shape


def test_fit_family_runs_on_reader_output(tmp_path):
    # end-to-end wiring: the reader's arrays feed fit_family without error.
    import axisym_family as af
    import fit_family as ff

    # synthesize a dump whose line spectra are broadband (not a family field);
    # we only check the pipeline executes and returns finite params.
    N, L = 4096, 4.0
    posf, xc, L = _uniform_line(N, L)
    rng = np.random.default_rng(2)
    # give v a red spectrum so the fit has a defined peak
    kfft = 2 * np.pi * np.fft.rfftfreq(N, d=L / N)
    amp = 1.0 / (1.0 + (kfft / 3.0) ** 2)
    phase = np.exp(2j * np.pi * rng.random(kfft.size))
    v = np.fft.irfft(amp * phase, n=N) * N
    u = 0.7 * v + 0.05 * rng.normal(0, 1, N)
    w = 0.7 * v + 0.05 * rng.normal(0, 1, N)
    odt_io.write_dump(os.path.join(tmp_path, "dmp_0000.dat"), posf, u, v, w, 1.0)

    k2, Eperp, E2 = odt_io.target_for_fit(str(tmp_path), strain=1.0, Nu=N,
                                          kmin=kfft[3], kmax=kfft[N // 8],
                                          nbin=40)
    assert k2.size <= 40
    res = ff.fit_family(k2, Eperp, E2)
    p = res.params
    assert np.isfinite([p.A0, p.c0, p.L2, p.Lperp]).all()
    assert p.A0 > 0 and p.L2 > 0 and p.Lperp > 0


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
