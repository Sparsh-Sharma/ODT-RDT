#!/usr/bin/env python3
"""Manuscript figure fig_level1a (main_v3.tex sec:level1a): verification
of the strain-coupled implementation in the rapid-distortion limit.
One two-panel canon figure: (a) component fractions — solver symbols on
the LRR moment trajectory, exact RDT for reference; (b) kinetic-energy
amplification.  Rerun 2026-09-06 on caro (case lvl1afig, deck
input/homogeneousStrain after fixing its committed merge-conflict
corruption); reference curves from fig_rdt_moments (verified).

Usage: python fig_level1a.py <dump_dir with dmp_*.dat>
Prints the verification numbers (max fraction deviation from LRR, kt
deviation) for the manuscript text.
"""
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
from figstyle_jfm import COL, FULL, panel, plt, save  # noqa: E402
from fig_rdt_moments import (EMAX, exact_rdt, integrate_model,  # noqa: E402
                             strain_tensor)


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
    posf = np.array(posf)
    u, v, w = np.array(u), np.array(v), np.array(w)
    faces = np.append(posf, -posf[0])   # centred domain
    dx = np.diff(faces)
    L = dx.sum()
    F = [q - (q * dx).sum() / L for q in (u, v, w)]
    R = np.array([[(F[a] * F[b] * dx).sum() / L for b in range(3)]
                  for a in range(3)])
    kt = 0.5 * np.trace(R)
    return t, np.array([R[i, i] / (2 * kt) for i in range(3)]), kt


def main():
    dump_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    rows = []
    for f in sorted(glob.glob(os.path.join(dump_dir, "dmp_*.dat"))):
        t, fr, kt = read_dump(f)
        if t is not None:
            rows.append((t, fr, kt))
    rows.sort(key=lambda r: r[0])
    e = np.array([r[0] for r in rows])
    fr = np.array([r[1] for r in rows])
    kt = np.array([r[2] for r in rows])
    print(f"read {len(rows)} dumps, e in [{e.min():.2f}, {e.max():.2f}]")

    A = strain_tensor("plane")
    te, Re = exact_rdt(A)
    tm, Rm = integrate_model(A, "LRR")
    ktm = 0.5 * (Rm[0, 0] + Rm[1, 1] + Rm[2, 2])
    fe = np.array([np.diag(Re[:, :, k]) for k in range(Re.shape[2])]).T
    fe /= fe.sum(0)
    fm = np.array([np.diag(Rm[:, :, k]) for k in range(Rm.shape[2])]).T
    fm /= fm.sum(0)

    # verification numbers for the text
    fm_at = np.array([[np.interp(ei, tm, fm[i]) for ei in e]
                      for i in range(3)])
    dev = np.max(np.abs(fm_at - fr.T))
    ktm_at = np.interp(e, tm, ktm / ktm[0])
    ktdev = np.max(np.abs(kt / kt[0] - ktm_at) / ktm_at)
    print(f"max |fraction - LRR| over all dumps/components: {dev:.2e}")
    print(f"max relative kt deviation from LRR: {ktdev:.2e}")

    CC = [COL["phi11"], COL["phi22"], COL["phi33"]]
    fig, (a, b) = plt.subplots(
        1, 2, figsize=(FULL, 2.3),
        gridspec_kw={"width_ratios": [1.35, 1.0]})
    for i in range(3):
        a.plot(te, fe[i], color=CC[i], lw=1.2)
        a.plot(tm, fm[i], color=CC[i], ls=(0, (6, 2)), lw=0.9)
        a.plot(e, fr[:, i], ls="none", marker="o", ms=2.8,
               mfc="none", mec=CC[i], mew=0.7)
    a.set_ylabel(r"$\overline{u_i^2}/2k_t$")
    from matplotlib.lines import Line2D
    hand = ([Line2D([0], [0], color=CC[i], lw=1.4) for i in range(3)]
            + [Line2D([0], [0], color="k", lw=1.2),
               Line2D([0], [0], color="k", ls=(0, (6, 2)), lw=0.9),
               Line2D([0], [0], color="k", ls="none", marker="o", ms=3,
                      mfc="none")])
    labs = [r"$\overline{u_1^2}$", r"$\overline{u_2^2}$",
            r"$\overline{u_3^2}$", "exact RDT", "LRR closure",
            "ODT solver"]
    a.legend(hand, labs, ncol=2, loc="upper left", fontsize=6.5,
             columnspacing=0.9, bbox_to_anchor=(0.01, 0.9))
    b.plot(tm, ktm / ktm[0], color="k", ls=(0, (6, 2)), lw=0.9,
           label="LRR closure")
    b.plot(e, kt / kt[0], ls="none", marker="o", ms=2.8, mfc="none",
           mec=COL["odt"], mew=0.7, label="ODT solver")
    b.set_ylabel(r"$k_t(e)/k_t(0)$")
    b.legend(loc="lower right")
    for j, ax in enumerate((a, b)):
        ax.set_xlabel(r"total strain $e = S\,t$")
        ax.set_xlim(0, EMAX)
        panel(ax, "ab"[j], x=0.02, y=0.15 if j == 0 else 0.97)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_level1a"))


if __name__ == "__main__":
    main()
