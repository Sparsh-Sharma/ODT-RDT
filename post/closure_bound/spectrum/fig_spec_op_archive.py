#!/usr/bin/env python3
"""Manuscript figure fig_spec_op (main_v3.tex sec:spec-op), rebuilt in
the JFM canon from the ARCHIVED vector figures
C_kv3em6_hires_op_{spectra,centroid}.pdf.

Provenance: the operating-point ensemble (nu=3e-6, 16000 cells, N=200,
V2 era) was produced in the pre-September workflow whose WSL environment
was lost; a caro reproduction stalled on the low-viscosity eddy-sampler
transient (deck fix committed: input/homogeneousStrainOP, Lmin pinned)
and is queued as a long-walltime task.  Meanwhile the plotted data are
recovered exactly from the archived vector PDFs.  Calibration is
validated internally: the archived rigid-translation reference in the
centroid panel must reproduce exp(e/2) (checked below), as must the
log-decade tick spacing of the spectra panel.
"""
import os
import sys

import fitz
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from figstyle_jfm import COL, FULL, panel, plt, save  # noqa: E402

FIGDIR = os.path.normpath(os.path.join(HERE, "..", "..", "..",
                                       "manuscript", "Figures"))
OLD = {(0.12, 0.47, 0.71): 0,   # blue was E1
       (0.84, 0.15, 0.16): 1,   # red  was E2
       (0.17, 0.63, 0.17): 2}   # green was E3
CC = [COL["phi11"], COL["phi22"], COL["phi33"]]


def polyline(D):
    pts = []
    for it in D["items"]:
        if it[0] == "l":
            if not pts:
                pts.append((it[1].x, it[1].y))
            pts.append((it[2].x, it[2].y))
    return np.array(pts)


def col_of(D):
    c = D.get("color") or D.get("fill")
    return tuple(np.round(c, 2)) if c else None


def logcal(words, axis):
    """Map pixel -> log10(value) from '10N' / '10-N' exponent labels."""
    px, ex = [], []
    for w in words:
        s = w[4].replace("−", "-")
        if s.startswith("10") and len(s) > 2:
            try:
                e = int(s[2:])
            except ValueError:
                continue
            cx = 0.5 * (w[0] + w[2])
            cy = 0.5 * (w[1] + w[3])
            if axis == "x" and 258 < w[1] < 272 and w[0] > 60:
                px.append(cx)
                ex.append(e)
            elif axis == "y" and w[0] < 50 and w[1] < 258:
                px.append(cy)
                ex.append(e)
    c = np.polyfit(px, ex, 1)
    span = np.ptp(np.polyval(c, px) - ex)
    assert span < 0.02, f"log calibration residual {span}"
    return lambda p: np.polyval(c, p)


def lincal(words, axis, xmax=50, ybot=260):
    px, v = [], []
    for w in words:
        try:
            val = float(w[4].replace("−", "-"))
        except ValueError:
            continue
        cx, cy = 0.5 * (w[0] + w[2]), 0.5 * (w[1] + w[3])
        if axis == "x" and w[1] > ybot:
            px.append(cx)
            v.append(val)
        elif axis == "y" and w[0] < xmax:
            px.append(cy)
            v.append(val)
    c = np.polyfit(px, v, 1)
    return lambda p: np.polyval(c, p)


def extract_spectra():
    p = fitz.open(os.path.join(FIGDIR,
                               "C_kv3em6_hires_op_spectra.pdf"))[0]
    words = p.get_text("words")
    fx = logcal(words, "x")
    fy = logcal(words, "y")
    curves, bands = {}, {}
    for D in p.get_drawings():
        c = col_of(D)
        if c not in OLD:
            continue
        bb = D["rect"]
        if bb.x1 - bb.x0 < 100:
            continue                       # legend handle
        pts = polyline(D)
        comp = OLD[c]
        rec = (10.0 ** fx(pts[:, 0]), 10.0 ** fy(pts[:, 1]))
        if D["type"] == "f":
            bands[comp] = rec
        else:
            curves[comp] = rec
    return curves, bands


def extract_centroid():
    p = fitz.open(os.path.join(FIGDIR,
                               "C_kv3em6_hires_op_centroid.pdf"))[0]
    words = p.get_text("words")
    fx = lincal(words, "x")
    fy = lincal(words, "y")
    out = {}
    for D in p.get_drawings():
        bb = D["rect"]
        c = col_of(D)
        if bb.x1 - bb.x0 < 100 or c in (None, (1.0, 1.0, 1.0)):
            continue
        pts = polyline(D)
        if pts.ndim != 2 or len(pts) < 8:
            continue
        rec = (fx(pts[:, 0]), fy(pts[:, 1]))
        if c == (0.0, 0.0, 0.0):
            out["rigid"] = rec
        elif D["type"] == "f":
            out["band"] = rec
        else:
            out["mean"] = rec
    # calibration check: archived rigid line == exp(e/2)
    x, y = out["rigid"]
    m = (x > 0.05) & (x < 3.95)
    err = np.max(np.abs(y[m] - np.exp(0.5 * x[m])))
    print(f"centroid calibration check |rigid - exp(e/2)|max = {err:.3f}")
    assert err < 0.12
    return out


def main():
    curves, bands = extract_spectra()
    cen = extract_centroid()

    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.3))
    for comp, lab in ((0, "$E_1$"), (1, "$E_2$"), (2, "$E_3$")):
        k, E = curves[comp]
        o = np.argsort(k)
        a.loglog(k[o], E[o], color=CC[comp], lw=0.8, label=lab)
        if comp in bands:
            kb, Eb = bands[comp]
            n = len(kb) // 2
            a.fill(kb, Eb, color=CC[comp], alpha=0.15, lw=0)
    kref = np.array([300.0, 6000.0])
    a.plot(kref, 1e-3 * (kref / kref[0]) ** (-5 / 3), "k:", lw=0.7)
    a.text(1500, 1.6e-3 * (1500 / kref[0]) ** (-5 / 3),
           r"$\kappa^{-5/3}$", fontsize=7)
    a.set_xlabel(r"wavenumber $\kappa_2$")
    a.set_ylabel(r"$E_i(\kappa_2)$")
    a.set_xlim(70, 1.3e5)
    a.set_ylim(1e-12, 3e-2)
    a.legend(loc="lower left")
    ee = np.linspace(0, 4, 100)
    b.plot(ee, np.exp(0.5 * ee), "k--", lw=0.9,
           label="rigid translation")
    x, y = cen["mean"]
    o = np.argsort(x)
    b.plot(x[o], y[o], color=COL["odt"], lw=1.1, label="ODT centroid")
    if "band" in cen:
        b.fill(cen["band"][0], cen["band"][1], color=COL["odt"],
               alpha=0.2, lw=0)
    b.set_xlabel(r"total strain $e$")
    b.set_ylabel(r"$\bar\kappa(e)/\bar\kappa(0)$")
    b.set_xlim(0, 4)
    b.set_ylim(0, 8)
    b.legend(loc="upper left")
    panel(a, "a", x=0.88)
    panel(b, "b", x=0.88, y=0.15)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_spec_op"))
    print("op centroid at e=4 (archived):",
          f"{np.interp(3.95, x[o], y[o]):.2f}")


if __name__ == "__main__":
    main()
