#!/usr/bin/env python3
"""Manuscript figures fig_LR_plane / fig_LR_axisym / fig_LR_b11_closure
(main_v3.tex sec:lr-validation), rebuilt in the JFM canon.

Data: ODT ensemble curves, standard-error bands and digitised L&R DNS
points are recovered from the archived vector figures (lr_extract.py —
exact to axis-calibration precision).  The analytic curves of the b11
figure (exact RDT, IP, LRR) are RECOMPUTED here from the Level-0
machinery under the L&R plane strain A = S*diag(0,-1,+1)/2 (S=1), and
the recomputed exact-RDT curve is cross-checked against the extracted
one — this validates both the calibration and the archived curve.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
from figstyle_jfm import COL, FULL, GOLD, GREEN, panel, plt, save  # noqa: E402
from fig_rdt_moments import (I3, exact_rdt, integrate_model,  # noqa: E402
                             production, rapid_IP, rapid_LRR)
from scipy.integrate import solve_ivp  # noqa: E402

D = np.load(os.path.join(HERE, "lr_extracted.npz"), allow_pickle=True)

OLD2COMP = {  # old figure colours -> component index
    (0.17, 0.63, 0.17): 0,   # green was b11
    (0.12, 0.47, 0.71): 1,   # blue was b22
    (0.84, 0.15, 0.16): 2,   # red was b33
}
CC = [COL["phi11"], COL["phi22"], COL["phi33"]]


def get(tag, kind, i, f):
    return D[f"{tag}_{kind}{i}_{f}"]


def col_of(tag, kind, i):
    return tuple(np.round(get(tag, kind, i, "color"), 2))


def comp_of(tag, kind, i):
    c = col_of(tag, kind, i)
    return OLD2COMP.get(c)


def band_edges(x, y):
    """Split a closed polygon into lower/upper edges on common x."""
    n = len(x) // 2
    return (x[:n], y[:n]), (x[n:][::-1], y[n:][::-1])


def lr_planeA():
    # L&R convention (V2 eq. total-strain): S = sqrt(S_ij S_ij / 2) = 1
    A = np.diag([0.0, -1.0, 1.0])
    S = 0.5 * (A + A.T)
    assert abs(np.sqrt(0.5 * np.sum(S * S)) - 1.0) < 1e-12
    return A


def b11_curves(emax=np.log(4.0), n=100):
    A = lr_planeA()
    te, Re = exact_rdt(A, emax=emax, nout=n)
    out = {}
    kt = 0.5 * (Re[0, 0] + Re[1, 1] + Re[2, 2])
    out["exact"] = (np.exp(te), Re[0, 0] / (2 * kt) - 1 / 3)
    for nm in ("IP", "LRR"):
        tm, Rm = integrate_model(A, nm, emax=emax, nout=n)
        ktm = 0.5 * (Rm[0, 0] + Rm[1, 1] + Rm[2, 2])
        out[nm] = (np.exp(tm), Rm[0, 0] / (2 * ktm) - 1 / 3)
    return out


def markers_by(tag, pred):
    xs, ys = [], []
    i = 0
    while f"{tag}_markers{i}_x" in D:
        if pred(i):
            xs.append(float(get(tag, "markers", i, "x")))
            ys.append(float(get(tag, "markers", i, "y")))
        i += 1
    o = np.argsort(xs)
    return np.array(xs)[o], np.array(ys)[o]


def fig_plane():
    fig, ax = plt.subplots(figsize=(0.8 * FULL, 2.6))
    names = ["$b_{11}$ (unstrained)", "$b_{22}$ (compressed, line)",
             "$b_{33}$ (stretched)"]
    mks = {0: "o", 1: "s", 2: "^"}
    # DNS b11 slow-fast range band (old green, 11-pt polygon)
    for i in range(4):
        if f"plane_bands{i}_x" not in D:
            continue
        x, y = get("plane", "bands", i, "x"), get("plane", "bands", i, "y")
        comp = comp_of("plane", "bands", i)
        if len(x) < 20:   # the DNS range band
            (xl, yl), (xu, yu) = band_edges(x, y)
            ax.fill_between(xl, yl, np.interp(xl, xu, yu),
                            color=CC[0], alpha=0.15, lw=0,
                            label=r"L&R DNS $b_{11}$ range "
                                  r"(slow--fast $S^{*}$)")
        else:             # ODT standard-error bands
            (xl, yl), (xu, yu) = band_edges(x, y)
            ax.fill_between(xl, yl, np.interp(xl, xu, yu),
                            color=CC[comp], alpha=0.25, lw=0)
    for i in range(3):
        comp = comp_of("plane", "curves", i)
        ax.plot(get("plane", "curves", i, "x"),
                get("plane", "curves", i, "y"),
                color=CC[comp], lw=1.2, label="ODT " + names[comp])
    for comp in range(3):
        x, y = markers_by(
            "plane", lambda i: comp_of("plane", "markers", i) == comp)
        ax.plot(x, y, ls="none", marker=mks[comp], ms=3.4, mfc="none",
                mec=CC[comp], mew=0.8,
                label=r"L&R DNS $b_{%d%d}$" % (comp + 1, comp + 1)
                      + (" (high $S^{*}$, RDT)" if comp == 0 else ""))
    ax.axhline(0, color="0.6", lw=0.5, ls=":")
    ax.set_xlabel(r"reference total strain $c=\exp\!\int S\,\mathrm{d}t$")
    ax.set_ylabel(r"$b_{ij}$")
    ax.legend(fontsize=6.3, loc="lower left", ncol=2,
              columnspacing=0.8)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_LR_plane"))


def fig_axisym():
    fig, ax = plt.subplots(figsize=(0.72 * FULL, 2.6))
    OLD = {(0.12, 0.47, 0.71): ("axis $b_{11}$", CC[0]),
           (0.84, 0.15, 0.16): (r"transverse $\frac{1}{2}(b_{22}+b_{33})$",
                                CC[1])}
    for i in range(2):
        c = col_of("axisym", "bands", i)
        lab, col = OLD[c]
        x, y = get("axisym", "bands", i, "x"), get("axisym", "bands", i,
                                                   "y")
        (xl, yl), (xu, yu) = band_edges(x, y)
        ax.fill_between(xl, yl, np.interp(xl, xu, yu), color=col,
                        alpha=0.25, lw=0,
                        label=("ODT $b_{22}$--$b_{33}$ spread"
                               if c == (0.84, 0.15, 0.16) else None))
    for i in range(2):
        lab, col = OLD[col_of("axisym", "curves", i)]
        ax.plot(get("axisym", "curves", i, "x"),
                get("axisym", "curves", i, "y"), color=col, lw=1.2,
                label="ODT " + lab)
    x, y = markers_by("axisym", lambda i: str(
        get("axisym", "markers", i, "shape")) == "re")
    ax.plot(x, y, "s", ms=3.4, mfc="none", mec=CC[0], mew=0.8,
            label=r"L&R DNS $b_{11}$")
    x, y = markers_by("axisym", lambda i: str(
        get("axisym", "markers", i, "shape")) == "l")
    ax.plot(x, y, "^", ms=3.4, mfc="none", mec=CC[1], mew=0.8,
            label=r"L&R DNS transverse")
    ax.axhline(0, color="0.6", lw=0.5, ls=":")
    ax.set_xlabel(r"reference total strain $c=\exp\!\int S\,\mathrm{d}t$")
    ax.set_ylabel(r"$b_{ij}$")
    ax.legend(fontsize=6.3, loc="lower left")
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_LR_axisym"))


def fig_b11():
    ana = b11_curves()
    # Calibration cross-check: the archived IP and LRR curves must match
    # the recomputation (they do, to <1e-3).  The archived "exact RDT"
    # curve is NOT used: it was constructed from the erroneous
    # kinematic-conservation argument (A_1k=0 => R_11 conserved), which
    # exact RDT violates through the rapid pressure-strain; the correct
    # exact curve, recomputed here, RISES to +0.114 at c=4 (discovered
    # 2026-09-07, propagates to the manuscript text).
    for i, nm in ((1, "IP"), (2, "LRR")):
        xg = get("b11", "curves", i, "x")
        o = np.argsort(xg)
        yr = np.interp(xg[o], *ana[nm])
        err = np.max(np.abs(yr - get("b11", "curves", i, "y")[o]))
        print(f"b11 calibration cross-check {nm}: "
              f"max |recomputed - extracted| = {err:.4f}")
    print("exact RDT b11 (recomputed) at c=1.5/2/3/4: "
          + " ".join(f"{np.interp(c, *ana['exact']):+.3f}"
                     for c in (1.5, 2, 3, 4)))

    fig, ax = plt.subplots(figsize=(0.72 * FULL, 2.6))
    # DNS range band (grey, small polygon)
    for i in range(2):
        x, y = get("b11", "bands", i, "x"), get("b11", "bands", i, "y")
        c = col_of("b11", "bands", i)
        (xl, yl), (xu, yu) = band_edges(x, y)
        if c == (0.75, 0.75, 0.75):
            ax.fill_between(xl, yl, np.interp(xl, xu, yu), color="0.8",
                            alpha=0.6, lw=0,
                            label=r"L&R DNS range (slow--fast $S^{*}$)")
        else:
            ax.fill_between(xl, yl, np.interp(xl, xu, yu),
                            color=COL["odt"], alpha=0.25, lw=0)
    ax.plot(*ana["exact"], color=COL["rdt"], ls=(0, (1, 1.5)), lw=1.3,
            label="exact linear RDT")
    ax.plot(*ana["IP"], color=GOLD, ls=(0, (4, 1.5, 1, 1.5)), lw=1.1,
            label=r"IP ($C_2=3/5$)")
    ax.plot(*ana["LRR"], color=GREEN, ls=(0, (6, 2)), lw=1.1,
            label="LRR-QI")
    ax.plot(get("b11", "curves", 3, "x"), get("b11", "curves", 3, "y"),
            color=COL["odt"], lw=1.3, label="present ODT ($N=1000$)")
    x, y = markers_by("b11", lambda i: True)
    ax.plot(x, y, "s", ms=3.4, mfc="none", mec="k", mew=0.8,
            label=r"L&R DNS (fast $S^{*}$)")
    ax.axhline(0, color="0.6", lw=0.5, ls=":")
    ax.set_xlabel(r"reference total strain $c=\exp\!\int S\,\mathrm{d}t$")
    ax.set_ylabel(r"$b_{11}$")
    ax.legend(fontsize=6.3, loc="upper left", ncol=2,
              columnspacing=0.8, frameon=True, framealpha=0.85,
              edgecolor="none")
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_LR_b11_closure"))


if __name__ == "__main__":
    fig_plane()
    fig_axisym()
    fig_b11()
