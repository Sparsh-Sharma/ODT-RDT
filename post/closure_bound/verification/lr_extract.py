#!/usr/bin/env python3
"""Extract the plotted data from the archived vector PDFs of the
Lee & Reynolds comparison figures (fig_LR_plane / fig_LR_axisym /
fig_LR_b11_closure).

Provenance: the ODT ensembles (N=1000) and the digitised L&R DNS points
were produced in the pre-September workflow whose WSL environment was
lost; the archived figures are vector graphics, so the plotted data are
recovered exactly (to axis-calibration precision, <0.5% of an axis
span).  The analytic curves (exact RDT, IP, LRR) are NOT taken from the
PDFs — fig_lr_paper.py recomputes them and uses the extracted exact-RDT
curve only as a cross-check of this calibration.

Writes lr_extracted.npz with, per figure, the calibrated polylines
(curves, band edges) and marker centres.
"""
import os

import fitz
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.normpath(os.path.join(HERE, "..", "..", "..",
                                       "manuscript", "Figures"))


def close(a, b, tol=1.5):
    return abs(a - b) <= tol


def calibrate(page):
    """Return (x_map, y_map) from tick marks + tick labels."""
    words = page.get_text("words")
    dr = page.get_drawings()
    # spine box: the longest black straight lines
    xs, ys = [], []
    for D in dr:
        if D["type"] not in ("s", "fs"):
            continue
        bb = D["rect"]
        if D.get("color") != (0.0, 0.0, 0.0):
            continue
        w, h = bb.x1 - bb.x0, bb.y1 - bb.y0
        if len(D["items"]) == 1 and D["items"][0][0] == "l":
            if w < 0.5 and 2 < h < 6:      # vertical tick
                xs.append((bb.x0, bb.y0, bb.y1))
            if h < 0.5 and 2 < w < 6:      # horizontal tick
                ys.append((bb.y0, bb.x0, bb.x1))
    # x ticks on the bottom: largest y0
    if not xs or not ys:
        raise RuntimeError("no ticks found")
    ymax = max(t[1] for t in xs)
    xticks = sorted({round(t[0], 1) for t in xs if t[1] > ymax - 8})
    xmin_ = min(t[1] for t in ys)
    yticks = sorted({round(t[0], 1) for t in ys if t[1] < xmin_ + 8})
    # numeric labels: x labels below the axes, y labels left of them
    # matplotlib may draw the minus sign as a PATH, so accept "0.15" etc.
    # and infer sign from position ordering afterwards.
    xlab, ylab = [], []
    for w in words:
        try:
            v = float(w[4].replace("−", "-"))
        except ValueError:
            continue
        cx, cy = 0.5 * (w[0] + w[2]), 0.5 * (w[1] + w[3])
        xlab.append((cx, cy, v))
        ylab.append((cx, cy, v))
    # associate: for each x tick, nearest label BELOW the axis
    ybot = max(xticks and [0] or [0])  # placeholder
    xpairs = []
    for xt in xticks:
        cand = [(abs(l[0] - xt), l) for l in xlab if l[1] > ymax]
        if cand:
            derr, lab = min(cand)
            if derr < 6:
                xpairs.append((xt, lab[2]))
    ypairs = []
    for yt in yticks:
        cand = [(abs(l[1] - yt), l) for l in ylab if l[0] < xmin_]
        if cand:
            derr, lab = min(cand)
            if derr < 6:
                ypairs.append((yt, lab[2]))
    xpairs.sort()
    ypairs.sort()
    # y labels lost their minus sign if drawn as a path: values must be
    # monotonically DECREASING with pixel y (data up = pixel down).
    vals = [v for _, v in ypairs]
    if len(set(vals)) < len(vals) or any(
            vals[i] < vals[i + 1] for i in range(len(vals) - 1)):
        # rebuild signs: assume uniform spacing, top tick has the largest
        # value; find the step from the two topmost distinct labels
        n = len(ypairs)
        step = abs(vals[0] - vals[1]) if n > 1 else 0.05
        top = max(vals)
        rebuilt = [(ypairs[i][0], top - i * step) for i in range(n)]
        # anchor: pixel row whose label reads 0.00 keeps value 0
        for py, v in ypairs:
            if v == 0.0:
                z = [r[0] for r in rebuilt]
                offs = 0.0 - rebuilt[int(np.argmin(np.abs(
                    np.array(z) - py)))][1]
                rebuilt = [(p, val + offs) for p, val in rebuilt]
                break
        ypairs = rebuilt
    px = np.array([p for p, _ in xpairs])
    vx = np.array([v for _, v in xpairs])
    py = np.array([p for p, _ in ypairs])
    vy = np.array([v for _, v in ypairs])
    cx = np.polyfit(px, vx, 1)
    cy = np.polyfit(py, vy, 1)
    return (lambda x: np.polyval(cx, x)), (lambda y: np.polyval(cy, y)), \
        (min(px), max(px))


def polyline(D):
    pts = []
    for it in D["items"]:
        if it[0] == "l":
            p1, p2 = it[1], it[2]
            if not pts:
                pts.append((p1.x, p1.y))
            pts.append((p2.x, p2.y))
    return np.array(pts)


def extract(pdfname):
    page = fitz.open(os.path.join(FIGDIR, pdfname))[0]
    fx, fy, (x0, x1) = calibrate(page)
    # legend handles: short (10-25 px) horizontal coloured line stubs;
    # any small marker within 30 px of one is a legend entry, not data
    handles = []
    for D in page.get_drawings():
        bb = D["rect"]
        w, h = bb.x1 - bb.x0, bb.y1 - bb.y0
        if (D["type"] in ("s", "fs") and 10 <= w <= 25 and h < 1.0
                and D.get("color") not in (None, (0.0, 0.0, 0.0))):
            handles.append((0.5 * (bb.x0 + bb.x1),
                            0.5 * (bb.y0 + bb.y1)))
    out = {"curves": [], "bands": [], "markers": []}
    for D in page.get_drawings():
        bb = D["rect"]
        cxp = 0.5 * (bb.x0 + bb.x1)
        col = D.get("color") or D.get("fill")
        if col == (1.0, 1.0, 1.0) or col is None:
            continue
        n = len(D["items"])
        wide = (bb.x1 - bb.x0) > 0.3 * (x1 - x0)
        small = (bb.x1 - bb.x0) <= 8 and (bb.y1 - bb.y0) <= 12
        if wide and n >= 8:
            pts = polyline(D)
            if len(pts) < 8:
                continue
            kind = "band" if (D["type"] in ("f", "fs")
                              and D.get("fill")) else "curve"
            # a matplotlib line saved as fs with fill==stroke and an
            # OPEN path is a curve; closed polygon (first==last) = band
            if kind == "band" and not (close(pts[0][0], pts[-1][0])
                                       and close(pts[0][1], pts[-1][1])):
                kind = "curve"
            rec = dict(color=col, dashes=D.get("dashes"),
                       x=fx(pts[:, 0]), y=fy(pts[:, 1]))
            out[kind + "s"].append(rec)
        elif small and cxp > x0 - 4 and D["type"] == "s":
            cyp = 0.5 * (bb.y0 + bb.y1)
            if any(abs(cxp - hx) < 30 and abs(cyp - hy) < 40
                   for hx, hy in handles):
                continue
            shape = D["items"][0][0]  # re / l / c
            out["markers"].append(dict(
                color=col, shape=shape,
                x=float(fx(cxp)), y=float(fy(cyp))))
    return out


def main():
    all_out = {}
    for tag, pdf in (("plane", "fig_LR_plane.pdf"),
                     ("axisym", "fig_LR_axisym.pdf"),
                     ("b11", "fig_LR_b11_closure.pdf")):
        res = extract(pdf)
        print(f"{tag}: {len(res['curves'])} curves, "
              f"{len(res['bands'])} bands, {len(res['markers'])} markers")
        for k, lst in res.items():
            for i, rec in enumerate(lst):
                for f, v in rec.items():
                    all_out[f"{tag}_{k}{i}_{f}"] = np.asarray(
                        v if v is not None else np.nan)
    np.savez(os.path.join(HERE, "lr_extracted.npz"), **all_out)
    print("wrote lr_extracted.npz")


if __name__ == "__main__":
    main()
