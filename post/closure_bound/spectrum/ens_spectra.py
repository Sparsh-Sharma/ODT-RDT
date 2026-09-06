#!/usr/bin/env python3
"""Ensemble spectral postprocessor for the spectrum-section campaigns
(sec:spec-bvc, sec:spec-op, sec:cmk-full in main_v3.tex).

For every realization of a case: read the dmp_*.dat series, resample to
a uniform grid, FFT the three components, and accumulate
  - component and total spectra at every dump time,
  - total-spectrum centroid and premultiplied-spectrum peak,
  - component (E1, E2) centroids,
  - kt(t), R11(t), R22(t) (for the energy budget), and the
    gradient-estimate dissipation 2*nu*sum_i <(du_i/dy)^2>.
Ensemble statistics: mean and geometric +-sigma (log-space) for spectra
and centroids; plain mean for energies.  Writes <case>_ens.npz.

Usage: python ens_spectra.py <data/CASE/data root> <out.npz> [Nu]
Runs on caro (heavy cases) or locally.
"""
import glob
import os
import re
import sys

import numpy as np

CENT_THRESH = 1e-3


def read_dump(fname):
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


def one_rlz(rdir, Nu):
    files = sorted(glob.glob(os.path.join(rdir, "dmp_*.dat")))
    out = []
    for f in files:
        t, posf, u, v, w = read_dump(f)
        if t is None:
            continue
        x0 = posf[0]
        L = -2.0 * x0
        faces = np.append(posf, -x0)
        xc = 0.5 * (faces[:-1] + faces[1:])
        dx = np.diff(faces)
        # moments (length-weighted)
        R = {}
        for nm, q in (("R11", u), ("R22", v), ("R33", w)):
            m = (q * dx).sum() / L
            R[nm] = ((q - m) ** 2 * dx).sum() / L
        kt = 0.5 * (R["R11"] + R["R22"] + R["R33"])
        # gradient dissipation (cell-centred differences)
        eps_g = 0.0
        for q in (u, v, w):
            g = np.gradient(q, xc)
            eps_g += 2.0 * ((g ** 2) * dx).sum() / L   # times nu later
        # spectra on uniform grid
        xu = x0 + (np.arange(Nu) + 0.5) * (L / Nu)
        k = 2.0 * np.pi * np.fft.rfftfreq(Nu, d=L / Nu)[1:]
        Es = []
        for q in (u, v, w):
            fu = np.interp(xu, xc, q, period=L)
            fu = fu - fu.mean()
            Es.append(2.0 * np.abs(np.fft.rfft(fu)[1:] / Nu) ** 2)
        Etot = Es[0] + Es[1] + Es[2]
        msk = Etot > CENT_THRESH * Etot.max()
        cent = (k[msk] * Etot[msk]).sum() / Etot[msk].sum()
        peak = k[np.argmax(k * Etot)]
        c1 = (k * Es[0]).sum() / Es[0].sum()
        c2 = (k * Es[1]).sum() / Es[1].sum()
        out.append(dict(t=t, L=L, kt=kt, R11=R["R11"], R22=R["R22"],
                        eps_g=eps_g, cent=cent, peak=peak, c1=c1, c2=c2,
                        k=k, E=Es, Etot=Etot))
    return out


def main():
    root, outnpz = sys.argv[1], sys.argv[2]
    Nu = int(sys.argv[3]) if len(sys.argv) > 3 else 4096
    rdirs = sorted(glob.glob(os.path.join(root, "data_*")))
    all_r = []
    for rd in rdirs:
        try:
            rows = one_rlz(rd, Nu)
        except Exception as ex:
            print(f"skip {rd}: {ex}")
            continue
        if rows:
            all_r.append(rows)
    n_dump = min(len(r) for r in all_r)
    N = len(all_r)
    print(f"{N} realizations x {n_dump} dumps")
    t = np.array([all_r[0][j]["t"] for j in range(n_dump)])

    def stack(key):
        return np.array([[r[j][key] for j in range(n_dump)]
                         for r in all_r])          # (N, n_dump)

    res = dict(t=t, N=N, L=stack("L")[0])
    for key in ("kt", "R11", "R22", "eps_g", "cent", "peak", "c1", "c2"):
        arr = stack(key)
        res[key + "_mean"] = arr.mean(0)
        res[key + "_sd"] = arr.std(0)
        if (arr > 0).all():
            lg = np.log(arr)
            res[key + "_gmean"] = np.exp(lg.mean(0))
            res[key + "_gsd"] = np.exp(lg.std(0))
    # ensemble spectra at first and last dump (common k grid per dump)
    for j, tag in ((0, "first"), (n_dump - 1, "last")):
        res[f"k_{tag}"] = all_r[0][j]["k"]
        for i in range(3):
            comp = np.array([r[j]["E"][i] for r in all_r])
            res[f"E{i + 1}_{tag}"] = comp.mean(0)
        tot = np.array([r[j]["Etot"] for r in all_r])
        res[f"Etot_{tag}"] = tot.mean(0)
        lg = np.log(np.maximum(tot, 1e-300))
        res[f"Etot_{tag}_gsd"] = np.exp(lg.std(0))
    np.savez(outnpz, **res)
    print("wrote", outnpz)


if __name__ == "__main__":
    main()
