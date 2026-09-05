#!/usr/bin/env python3
"""Inspect the gateA_S1 pilot realization: spectra at each dump + a trial
vK-family fit at e=0 (post-precursor) and e=2.

    python3 pilot_gateA_look.py <dump_dir>
"""
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import odt_io                                          # noqa: E402
import fit_family as ff                                # noqa: E402

root = sys.argv[1]
files = sorted(glob.glob(os.path.join(root, "dmp_*.dat")))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
cols = plt.cm.viridis(np.linspace(0, 0.9, len(files)))
fits = {}
for fn, col in zip(files, cols):
    t, posf, u, v, w = odt_io.read_dump(fn)
    k2, E1, E2, E3 = odt_io.component_spectra(posf, u, v, w, Nu=8192)
    e = max(t - 0.4, 0.0)
    ax1.loglog(k2, E2, color=col, lw=0.9, label=f"e={e:.1f} (n={posf.size})")
    if e in (0.0, 2.0):
        kmin = 3.0 * k2[0]
        sel = (k2 >= kmin) & (k2 <= 3000.0)
        k2b, (Epb, E2b) = odt_io.log_bin(k2[sel],
                                         [0.5*(E1[sel]+E3[sel]), E2[sel]], 36)
        res = ff.fit_family(k2b, Epb, E2b)
        p = res.params
        fits[e] = (res, k2b, E2b)
        E1m, E2m = ff.line_spectra(p, k2b)
        ax2.loglog(k2b, E2b, "o", ms=3, color=col, label=f"e={e:.1f} data")
        ax2.loglog(k2b, E2m, "-", color=col,
                   label=f"fit: L2={p.L2:.0f} Lp={p.Lperp:.0f} "
                         f"c0={p.c0:.2f} cost={res.cost:.1f}")
        print(f"e={e}: A0={p.A0:.3e} c0={p.c0:.3f} L2={p.L2:.2f} "
              f"Lperp={p.Lperp:.2f} ratio={p.Lperp/p.L2:.2f} "
              f"cost={res.cost:.2f} nbin={k2b.size}")
ax1.axvline(2*np.pi/0.0002/2, color="k", lw=0.6, ls=":")
ax1.set_title("pilot $E_2(k_2)$, all dumps")
ax2.set_title("single-realization trial fits")
for ax in (ax1, ax2):
    ax.set_xlabel("$k_2$")
    ax.legend(fontsize=6)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_gateA_pilot.png"), dpi=180)
print("saved fig_gateA_pilot.png")
