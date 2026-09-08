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
from figstyle_jfm import FULL, panel, plt, save  # noqa: E402
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

    # monochrome: source -> style (solid exact / dashed LRR / open
    # circles solver, drawn sparse so the LRR dash shows through);
    # component identified by direct text labels, not the legend.
    fig, (a, b) = plt.subplots(
        1, 2, figsize=(FULL, 2.3),
        gridspec_kw={"width_ratios": [1.35, 1.0]})
    # monochrome: component -> marker on the model band (o/s/^, open =
    # ODT solver, lying on the dashed LRR line); source -> line style
    # (solid exact / dashed LRR).  Exact solids are text-labelled so the
    # crossing u2/u3 lines stay identifiable.
    # monochrome: shape -> component (o/s/^ = u1/u2/u3); fill -> source
    # (filled = exact RDT on the solid line, open = ODT solver on the
    # dashed LRR line).  Crossing-proof, no in-plot text needed.
    MK = ["o", "s", "^"]                       # u_1, u_2, u_3
    for i in range(3):
        a.plot(te, fe[i], color="k", lw=1.1, marker=MK[i],
               markevery=(8 * i + 6, 40), ms=3.0, mfc="k", mec="k")
        a.plot(tm, fm[i], color="k", ls=(0, (6, 2)), lw=1.0)
        a.plot(e, fr[:, i], ls="none", marker=MK[i], markevery=2,
               ms=3.4, mfc="none", mec="k", mew=0.8)
    a.set_ylabel(r"$\overline{u_i^2}/2k_t$")
    b.plot(tm, ktm / ktm[0], color="k", ls=(0, (6, 2)), lw=1.0,
           label="LRR closure")
    b.plot(e, kt / kt[0], ls="none", marker="o", markevery=2, ms=3.2,
           mfc="none", mec="k", mew=0.8, label="SC-ODT (solver)")
    b.set_ylabel(r"$k_t(e)/k_t(0)$")
    for j, ax in enumerate((a, b)):
        ax.set_xlabel(r"total strain $e = S\,t$")
        ax.set_xlim(0, EMAX)
        panel(ax, "ab"[j], x=0.03, y=0.96)
    from matplotlib.lines import Line2D
    # panel (a): compact SOURCE legend (far right, clear band where the
    # u_1 solver curve has dropped to ~0.1); components identified by
    # direct labels on the three well-separated solver bands.
    hand = [Line2D([0], [0], color="k", lw=1.1, marker="o", mfc="k",
                   ms=3.5),
            Line2D([0], [0], color="k", ls=(0, (6, 2)), lw=1.0),
            Line2D([0], [0], color="k", ls="none", marker="o",
                   mfc="none", mew=0.9, ms=4)]
    la = a.legend(hand, ["exact RDT", "LRR closure", "SC-ODT (solver)"],
                  loc="center right", bbox_to_anchor=(0.99, 0.35),
                  fontsize=6.2, handlelength=1.7, labelspacing=0.3,
                  borderpad=0.4)
    la.get_frame().set_edgecolor("0.7")
    for lab, xv, yv in ((r"$\overline{u_2^2}$", 0.85, 0.505),
                        (r"$\overline{u_3^2}$", 1.15, 0.285),
                        (r"$\overline{u_1^2}$", 0.85, 0.135)):
        a.text(xv, yv, lab, fontsize=8, ha="center", va="center")
    b.legend(loc="upper left", bbox_to_anchor=(0.05, 0.86),
             fontsize=6.8, borderpad=0.4)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_level1a"))


if __name__ == "__main__":
    main()
