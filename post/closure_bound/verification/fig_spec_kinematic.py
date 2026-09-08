#!/usr/bin/env python3
"""Manuscript figure fig_spec_kinematic (main_v3.tex sec:spec-rdt):
compression without cascade — eddy events suppressed, inviscid, domain
dilatation on.  (a) component line spectra at increasing total strain
(lighter to darker); (b) spectral-centroid migration against the
kinematic prediction exp(-A22 e).  Rerun 2026-09-06 on caro (case
hsA2fig, deck input/homogeneousStrainA2); replaces the V2-era
fig_spectrum_hs__spectra/__centroid pair (generator lost with the old
WSL distro).

Usage: python fig_spec_kinematic.py <case dump dir>
Prints the centroid ratios at the final strain for the text.
"""
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from figstyle_jfm import FULL, panel, plt, save  # noqa: E402

A22 = -0.5
NUNIFORM = 4096
STRAINS = [0.0, 1.0, 2.0, 3.0, 3.9]
CENT_THRESH = 1e-3
# monochrome: component -> line style, strain -> grey ramp (light -> black)
LS = ["-", (0, (6, 2)), (0, (1, 1.4))]          # E_1, E_2, E_3
MK = ["o", "s", "^"]


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


def spectra_of(posf, fields, Nu=NUNIFORM):
    x0 = posf[0]
    L = -2.0 * x0
    faces = np.append(posf, -x0)
    xc = 0.5 * (faces[:-1] + faces[1:])
    xu = x0 + (np.arange(Nu) + 0.5) * (L / Nu)
    k = 2.0 * np.pi * np.fft.rfftfreq(Nu, d=L / Nu)
    out = []
    for f in fields:
        fu = np.interp(xu, xc, f, period=L)
        fu = fu - fu.mean()
        out.append(2.0 * np.abs(np.fft.rfft(fu) / Nu) ** 2)
    return k[1:], [E[1:] for E in out]


def centroid(k, E):
    m = E > CENT_THRESH * E.max()
    return (k[m] * E[m]).sum() / E[m].sum()


def main():
    case = sys.argv[1] if len(sys.argv) > 1 else "."
    rows = []
    for f in sorted(glob.glob(os.path.join(case, "dmp_*.dat"))):
        t, posf, u, v, w = read_dump(f)
        if t is not None:
            rows.append((t, posf, u, v, w))
    rows.sort(key=lambda r: r[0])
    times = np.array([r[0] for r in rows])
    L0 = -2.0 * rows[0][1][0]
    Lf = -2.0 * rows[-1][1][0]
    print(f"{len(rows)} dumps, e in [{times[0]:.2f},{times[-1]:.2f}], "
          f"L/L0 = {Lf / L0:.4f} (exact {np.exp(A22 * times[-1]):.4f})")

    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.3))
    epts, cent3 = [], {0: [], 1: [], 2: []}
    for e_t in STRAINS:
        j = int(np.argmin(np.abs(times - e_t)))
        _, posf, u, v, w = rows[j]
        k, Es = spectra_of(posf, [u, v, w])
        grey = 0.62 * (1.0 - e_t / max(STRAINS))   # 0.62 (e=0) -> 0 (black)
        for i in range(3):
            a.loglog(k, Es[i], color=str(grey), ls=LS[i], lw=0.85)
        epts.append(times[j])
        for i in range(3):
            cent3[i].append(centroid(k, Es[i]))
    epts = np.array(epts)
    for i in range(3):
        r = np.array(cent3[i]) / cent3[i][0]
        print(f"centroid ratio comp {i + 1} at e={epts[-1]:.1f}: "
              f"{r[-1]:.2f} (exact {np.exp(-A22 * epts[-1]):.2f})")
    a.set_xlabel(r"wavenumber $\kappa_2$")
    a.set_ylabel(r"$E_i(\kappa_2)$")
    a.set_xlim(6, 3e3)
    a.set_ylim(1e-5, 2)
    from matplotlib.lines import Line2D
    a.legend([Line2D([0], [0], color="k", ls=LS[i], lw=1.3)
              for i in range(3)],
             [r"$E_1$", r"$E_2$", r"$E_3$"], loc="upper left",
             bbox_to_anchor=(0.0, 1.0),
             title=r"light $\to$ black: $e\!=\!0\to3.9$")
    a.get_legend().get_title().set_fontsize(6.5)
    ee = np.linspace(0, max(STRAINS), 100)
    b.plot(ee, np.exp(-A22 * ee), "k--", lw=0.9,
           label=r"kinematic $\mathrm{e}^{-A_{22}e}$")
    for i in range(3):
        b.plot(epts, np.array(cent3[i]) / cent3[i][0], color="k",
               ls="-", lw=0.7, marker=MK[i], ms=3.6, mfc="none", mew=0.8,
               label=r"$\bar\kappa_%d$" % (i + 1))
    b.set_xlabel(r"total strain $e = S\,t$")
    b.set_ylabel(r"$\bar\kappa(e)/\bar\kappa(0)$")
    b.legend(loc="upper left")
    panel(a, "a", x=0.88)
    panel(b, "b", x=0.9, y=0.15)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_spec_kinematic"))


if __name__ == "__main__":
    main()
