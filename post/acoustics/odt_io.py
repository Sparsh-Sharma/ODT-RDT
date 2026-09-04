#!/usr/bin/env python3
r"""
Read homogeneous-strain ODT dumps and reduce them to the component line spectra
E_i(k2) that fit_family.fit_family targets (paper 2, Gate A, Session 2).

The reduction matches spectrum_diagnostic2.py exactly (interpolate the adapted,
dilated line onto a uniform grid over the actual domain L(e) = -2*posf[0], remove
the mean, rFFT each velocity component), with ONE convention change for fitting:
the returned spectra are spectral DENSITIES E_i(k2) [= |FFT|^2 * 2 / dk], so that

    int E_i(k2) dk2  ~=  <u_i^2> = R_ii,

the same normalization axisym_family.E_component / R_ii use.  (The raw per-mode
PSD used in spectrum_diagnostic2.py differs from this by the constant dk = 2*pi/L,
which a log-space fit would otherwise absorb into A0; converting here keeps the
fitted A0 -- hence k_t -- physically meaningful.)

Dump columns (whitespace-separated, '#'-comment header carrying `time = ...`):
    idx   pos_face   u(=u1)   v(=u2)   w(=u3)
so component 2 (the FFT of v) is the upwash line spectrum E_2(k2).

DATA-INDEPENDENT MODULE: it reads whatever path you give it; the repo ships no
homogeneous-strain dumps, so tests here drive it with synthetic dumps written by
`write_dump`.
"""

import glob
import os
import re

import numpy as np


# ----------------------------------------------------------------------
# Low-level dump IO
# ----------------------------------------------------------------------
def read_dump(fname):
    """Read one dump.  Returns (t, posf, u, v, w) with posf the cell faces."""
    t = None
    posf, u, v, w = [], [], [], []
    with open(fname) as f:
        for line in f:
            if line.startswith("#"):
                m = re.search(r"time\s*=\s*([-\d.eE+]+)", line)
                if m:
                    t = float(m.group(1))
                continue
            if not line.strip():
                continue
            c = line.split()
            posf.append(float(c[1]))
            u.append(float(c[2]))
            v.append(float(c[3]))
            w.append(float(c[4]))
    return t, np.array(posf), np.array(u), np.array(v), np.array(w)


def write_dump(fname, posf, u, v, w, t=0.0):
    """Write a dump in the read_dump format (used for tests / synthetic data)."""
    posf = np.asarray(posf, float)
    with open(fname, "w") as f:
        f.write(f"# time = {t:.8e}\n")
        f.write("# idx  pos  u  v  w\n")
        for i in range(posf.size):
            f.write(f"{i:6d}  {posf[i]:.8e}  {u[i]:.8e}  "
                    f"{v[i]:.8e}  {w[i]:.8e}\n")


def find_dumps(root):
    """Locate dump files under `root` (handles the run/.../data_*/dmp_*.dat and
    flat layouts, as in spectrum_diagnostic2.py)."""
    if os.path.isfile(root):
        return [root]
    for cand in (os.path.join(root, "data"), root):
        fs = (glob.glob(os.path.join(cand, "data_*", "dmp_*.dat"))
              or glob.glob(os.path.join(cand, "dmp_*.dat")))
        if fs:
            return sorted(fs)
    return sorted(glob.glob(os.path.join(root, "*.dat")))


# ----------------------------------------------------------------------
# Line-FFT reduction to component spectral densities E_i(k2)
# ----------------------------------------------------------------------
def component_spectra(posf, u, v, w, Nu=2048):
    """Reduce one line to component spectral DENSITIES.

    Returns (k2, E1, E2, E3) with k2 the positive wavenumbers along the ODT
    axis (k=0 dropped) and E_i spectral densities (int E_i dk2 ~= <u_i^2>).
    """
    x0 = posf[0]
    L = -2.0 * x0
    faces = np.append(posf, -x0)
    xc = 0.5 * (faces[:-1] + faces[1:])
    xu = x0 + (np.arange(Nu) + 0.5) * (L / Nu)
    k = 2.0 * np.pi * np.fft.rfftfreq(Nu, d=L / Nu)
    dk = 2.0 * np.pi / L  # uniform mode spacing
    out = []
    for f in (u, v, w):
        fu = np.interp(xu, xc, f, period=L)
        fu = fu - fu.mean()
        fh = np.fft.rfft(fu) / Nu
        psd = 2.0 * np.abs(fh) ** 2      # one-sided per-mode power
        out.append(psd / dk)             # -> spectral density
    return k[1:], out[0][1:], out[1][1:], out[2][1:]


# ----------------------------------------------------------------------
# Fit-target loader
# ----------------------------------------------------------------------
def load_fit_target(root, strain=None, smag=1.0, Nu=2048):
    """Load component line spectra to fit, at (or nearest) a target strain.

    Parameters
    ----------
    root   : a dump file, or a directory searched by find_dumps.
    strain : target total strain e = smag * t.  If None, the last dump is used.
    smag   : strain magnitude S so that e = smag * t (cf. spectrum_diagnostic2).
    Nu     : uniform resampling resolution for the FFT.

    Returns
    -------
    dict with keys k2, E1, E2, E3 (spectral densities), plus e, t, L, file.
    Note E1 and E3 are the perpendicular components (equal in the model);
    E2 is the upwash line spectrum.
    """
    files = find_dumps(root)
    if not files:
        raise FileNotFoundError(f"no dumps found under {root!r}")

    recs = []
    for fn in files:
        t, posf, u, v, w = read_dump(fn)
        if t is None or posf.size == 0:
            continue
        recs.append((t, fn, posf, u, v, w))
    if not recs:
        raise ValueError(f"no readable dumps under {root!r}")
    recs.sort(key=lambda r: r[0])

    if strain is None:
        pick = recs[-1]
    else:
        pick = min(recs, key=lambda r: abs(smag * r[0] - strain))
    t, fn, posf, u, v, w = pick
    k2, E1, E2, E3 = component_spectra(posf, u, v, w, Nu=Nu)
    return {
        "k2": k2, "E1": E1, "E2": E2, "E3": E3,
        "e": smag * t, "t": t, "L": -2.0 * posf[0], "file": fn,
    }


def log_bin(k, values, nbin):
    """Average `values` (list or array) onto `nbin` log-spaced k bins.

    Returns (k_binned, [values_binned...]).  Empty bins are dropped.  This
    denoises the dense FFT spectrum and gives the fitter a modest, evenly
    log-spread set of points.
    """
    k = np.asarray(k, float)
    single = np.ndim(values) == 1 and not isinstance(values, (list, tuple))
    arrs = [np.asarray(values, float)] if single else [np.asarray(v, float)
                                                       for v in values]
    edges = np.geomspace(k.min(), k.max(), nbin + 1)
    idx = np.clip(np.digitize(k, edges) - 1, 0, nbin - 1)
    kb, out = [], [[] for _ in arrs]
    for b in range(nbin):
        m = idx == b
        if not np.any(m):
            continue
        kb.append(np.exp(np.mean(np.log(k[m]))))
        for j, a in enumerate(arrs):
            out[j].append(np.mean(a[m]))
    kb = np.array(kb)
    outs = [np.array(o) for o in out]
    return (kb, outs[0]) if single else (kb, outs)


def ensemble_fit_target(root, strain=None, smag=1.0, Nu=2048):
    """Like load_fit_target, but averages the component spectra over ALL
    realization directories (data_00000, data_00001, ...) under `root`,
    picking in each the dump nearest the target strain.

    The dilatation history is deterministic, so all realizations share the
    same domain length and hence the same FFT wavenumber grid at a given
    time; spectra are averaged pointwise.  Returns the load_fit_target dict
    plus 'nrlz'.
    """
    files = find_dumps(root)
    if not files:
        raise FileNotFoundError(f"no dumps found under {root!r}")
    groups = {}
    for fn in files:
        groups.setdefault(os.path.dirname(fn), []).append(fn)

    acc, k2ref, meta, nrlz = None, None, None, 0
    for _, fns in sorted(groups.items()):
        recs = []
        for fn in fns:
            t, posf, u, v, w = read_dump(fn)
            if t is None or posf.size == 0:
                continue
            recs.append((t, fn, posf, u, v, w))
        if not recs:
            continue
        recs.sort(key=lambda r: r[0])
        pick = recs[-1] if strain is None else min(
            recs, key=lambda r: abs(smag * r[0] - strain))
        t, fn, posf, u, v, w = pick
        k2, E1, E2, E3 = component_spectra(posf, u, v, w, Nu=Nu)
        if acc is None:
            k2ref = k2
            acc = [np.zeros_like(E1) for _ in range(3)]
            meta = {"e": smag * t, "t": t, "L": -2.0 * posf[0], "file": fn}
        elif not np.allclose(k2, k2ref, rtol=1e-6):
            # differing domain length (shouldn't happen): interpolate in log k
            E1, E2, E3 = (np.interp(np.log(k2ref), np.log(k2), E)
                          for E in (E1, E2, E3))
        for a, E in zip(acc, (E1, E2, E3)):
            a += E
        nrlz += 1
    out = {"k2": k2ref, "E1": acc[0] / nrlz, "E2": acc[1] / nrlz,
           "E3": acc[2] / nrlz, "nrlz": nrlz}
    out.update(meta)
    return out


def target_for_fit(root, strain=None, smag=1.0, Nu=2048,
                   kmin=None, kmax=None, nbin=None):
    """Convenience wrapper returning (k2, E_perp, E2) ready for
    fit_family.fit_family, with E_perp = (E1+E3)/2 and an optional [kmin,kmax]
    band to exclude the low-k finite-domain artefact and the high-k FFT floor.
    If `nbin` is given, the band is log-binned to nbin points (denoise + speed).
    """
    d = load_fit_target(root, strain=strain, smag=smag, Nu=Nu)
    k2 = d["k2"]
    Eperp = 0.5 * (d["E1"] + d["E3"])
    E2 = d["E2"]
    sel = np.ones_like(k2, dtype=bool)
    if kmin is not None:
        sel &= k2 >= kmin
    if kmax is not None:
        sel &= k2 <= kmax
    k2, Eperp, E2 = k2[sel], Eperp[sel], E2[sel]
    if nbin is not None:
        k2, (Eperp, E2) = log_bin(k2, [Eperp, E2], nbin)
    return k2, Eperp, E2
