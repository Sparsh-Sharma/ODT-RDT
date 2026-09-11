#!/usr/bin/env python3
"""CANDIDATE money-shot figure for Part III (Sparsh 2026-09-16): the cost
of the frozen-inflow assumption in one panel.

Signed far-field SPL prediction error against total distortion strain e,
judged by the best available reference at each strain:

  laboratory window (e_eff ~ 0.03-0.16): the eight published experiments
    of fig_expval (biases from metrics_expval.txt, strains from
    cases.py); flat plates filled (clean), NACA0012 open (its common
    overprediction is the finite-thickness effect, identical for every
    inflow).
  beyond it (e = 0.5-2): the exact-RDT-distorted spectrum of the model's
    relaxed state (the reference DNS validates at rapid rates); band-mean
    error over 100-1000 Hz at the fig_amiet_rapid configuration, model
    bands = +-10 log10(1+rms) where the measured rms distance exists
    (e = 1, 2).

Series: frozen von Karman (standard practice), computed SC-ODT and
clock-relaxed ODT inflows.  Monochrome, JFM canon."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import fig_amiet_rapid as far  # noqa: E402
from figstyle_jfm import FULL, plt, save  # noqa: E402

# ---- laboratory window: biases [dB] from metrics_expval.txt, strains
#      (e_ctr, e_lo, e_hi) from post/acoustics/expval/cases.py
PLATES = [  # (case, e_ctr, e_lo, e_hi, vk, odt_std, odt_fix)
    ("NA15_60", 0.0363, 0.022, 0.056, -1.74, -1.86, -1.86),
    ("BA22_19", 0.0477, 0.028, 0.078, +1.52, +1.33, +1.32),
    ("BA22_27", 0.0477, 0.028, 0.078, -0.24, -0.40, -0.40),
    ("BA22_32", 0.0477, 0.028, 0.078, +0.28, +0.12, +0.11),
    ("BA19_32", 0.0477, 0.028, 0.078, +0.08, -0.08, -0.09),
]
NACA = [
    ("PA40", 0.1602, 0.088, 0.273, +5.70, +5.35, +5.33),
    ("PA60", 0.1606, 0.088, 0.274, +1.80, +1.42, +1.39),
    ("PA90", 0.1637, 0.090, 0.278, +5.12, +4.70, +4.66),
]

ES_APP = [0.5, 1.0, 2.0]
FBAND = (100.0, 1000.0)

# per-series x jitter (multiplicative, log axis) so coincident lab points
# stay legible
JIT = {"vk": 1.0, "std": 1.055, "fix": 1.115}
STYLE = {  # colour, marker, line style, label
    "vk": ("k", "s", (0, (5, 2)), "frozen von Kármán"),
    "std": ("0.45", "o", "-", "computed, SC-ODT"),
    "fix": ("0.45", "^", (0, (1, 1.4)), "computed, clock-relaxed ODT"),
}


def app_errors():
    """Band-mean SPL error vs the exact-RDT reference at ES_APP."""
    fits_s = far.read_fits("S20")
    fits_f = far.read_fits("S20_RCS1")
    sp_s = np.load(os.path.join(far.HERE, "spectra_gateA_S20.npz"),
                   allow_pickle=True)
    sp_f = np.load(os.path.join(far.HERE, "spectra_gateA_S20_RCS1.npz"),
                   allow_pickle=True)
    ke0_s, ke0_f = fits_s[0.0][0], fits_f[0.0][0]
    ke_phys = far.CVK / far.LAMBDA
    w2 = far.UP ** 2
    f = np.geomspace(30.0, 2000.0, 40)
    m = (f >= FBAND[0]) & (f <= FBAND[1])
    out = {k: [] for k in ("vk", "std", "fix")}
    for e in ES_APP:
        g22 = far.variance_m(e, ke0_s) / far.variance_m(0.0, ke0_s)
        spl_t = far.spl_curve(f, e, ke0_s * np.exp(e / 2.0),
                              ke0_s / ke_phys, g22, w2)
        spl_v = far.spl_curve(f, 0.0, ke0_s, ke0_s / ke_phys, 1.0, w2)
        out["vk"].append(np.mean((spl_v - spl_t)[m]))
        for tag, fits, sp, ke0 in (("std", fits_s, sp_s, ke0_s),
                                   ("fix", fits_f, sp_f, ke0_f)):
            lev = far.r22(sp, e) / far.r22(sp, 0.0)
            spl_m = far.spl_curve(f, e, fits[e][0], ke0 / ke_phys, lev, w2)
            out[tag].append(np.mean((spl_m - spl_t)[m]))
    return out


def main():
    app = app_errors()
    fig, ax = plt.subplots(figsize=(0.8 * FULL, 2.75))

    # regions: laboratory window vs beyond it
    ax.axvspan(0.020, 0.29, color="0.94", zorder=0)
    ax.axhline(0.0, color="0.75", lw=0.6, zorder=1)
    ax.text(0.076, 6.55, "laboratory window\n(reference: measurement)",
            fontsize=7, ha="center", va="top", color="0.25")
    ax.text(1.05, 6.55, "beyond the laboratory window\n"
            "(reference: exact RDT, DNS-anchored)",
            fontsize=7, ha="center", va="top", color="0.25")

    for key in ("vk", "std", "fix"):
        col, mk, ls, lab = STYLE[key]
        j = JIT[key]
        # flat plates: filled; NACA0012: open (thickness-contaminated)
        col_i = {"vk": 4, "std": 5, "fix": 6}[key]
        for rows, filled in ((PLATES, True), (NACA, False)):
            x = np.array([r[1] for r in rows]) * j
            xe = np.array([[r[1] - r[2] for r in rows],
                           [r[3] - r[1] for r in rows]]) * j
            yv = [r[col_i] for r in rows]
            ax.errorbar(x, yv, xerr=xe, fmt=mk, color=col, ms=3.2,
                        mfc=(col if filled else "white"), mec=col,
                        ecolor="0.82", elinewidth=0.5, capsize=0,
                        ls="none", zorder=4)
        # application strains: connected curve (+ measured-rms band)
        ax.plot(ES_APP, app[key], color=col, ls=ls, marker=mk, ms=3.2,
                mfc=col, lw=1.1, zorder=5, label=lab)
        if key in far.RMS:
            eb = [e for e in ES_APP if e in far.RMS[key]]
            db = np.array([10 * np.log10(1 + far.RMS[key][e]) for e in eb])
            yb = np.array([app[key][ES_APP.index(e)] for e in eb])
            ax.fill_between(eb, yb - db, yb + db, color="0.85",
                            alpha=0.55, lw=0, zorder=2)

    # the finite-thickness cluster, called out once
    ax.annotate("NACA0012: finite-thickness\noverprediction, common\n"
                "to every inflow (open)",
                xy=(0.150, 5.0), xytext=(0.036, 3.05), fontsize=7,
                color="0.25", ha="left", va="top",
                arrowprops=dict(arrowstyle="-", color="0.5", lw=0.5,
                                shrinkB=4))
    ax.annotate("the omitted term:\ncost of freezing the inflow",
                xy=(1.90, 3.15), xytext=(0.47, 4.35), fontsize=7,
                ha="center",
                arrowprops=dict(arrowstyle="-|>", color="k", lw=0.7,
                                shrinkB=3))

    ax.set_xscale("log")
    ax.set_xlim(0.020, 3.2)
    ax.set_ylim(-4.4, 6.9)
    ax.set_xticks([0.03, 0.1, 0.3, 1, 2])
    ax.set_xticklabels(["0.03", "0.1", "0.3", "1", "2"])
    ax.set_xlabel(r"total distortion strain $e$")
    ax.set_ylabel(r"prediction error [dB]")
    ax.legend(loc="lower left", fontsize=7)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_frozen_cost"))


if __name__ == "__main__":
    main()
