#!/usr/bin/env python3
"""Digitise Zusi & Perot 2013 (Phys. Fluids 25, 110819), FIG. 10(a), p. 12:
diagonal anisotropy components b_ij of IC3 under uniform plane strain,
against S*(t - t0) = total strain e, for several strain rates.  Vector
figure: every curve is a PDF path; points are the segment end-points.

Straining lasts until e = 0.5 in every run (S*Delta t = 0.5, their Sec.
II); beyond that the curves are the return to isotropy in the same
S-scaled time, so only e <= 0.5 is a strained state.  We keep
e in [0, EMAX] and write one CSV per (component, path).

Component convention (theirs -> ours, from the signs in their figure and
the RDT physics: compressed direction gains energy):
  their b22 (green, positive)  = compressed  -> our b22
  their b11 (red,   negative)  = stretched   -> our b33
  their b33 (blue,  ~0)        = neutral     -> our b11

Calibration: linear fits pdf->data on the numeric tick labels of panel
(a); the fit residual is asserted below TOL so a mis-read tick cannot
pass silently.  Usage: python extract_zp2013.py  (writes zp2013_*.csv)
"""
import csv
import os

import fitz
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "ZusiPerot2013.pdf")
PAGE = 12
EMAX = 0.55
TOL = 0.015                       # calibration residual, data units
COMP = {(0.05, 0.51, 0.25): ("b22", "compressed -> our b22"),
        (0.93, 0.13, 0.14): ("b11", "stretched  -> our b33"),
        (0.23, 0.33, 0.64): ("b33", "neutral    -> our b11")}


def calib(words):
    xs, xv, ys, yv = [], [], [], []
    for w in words:
        s = w[4].replace("−", "-").replace("–", "-")
        try:
            v = float(s)
        except ValueError:
            continue
        cx, cy = 0.5 * (w[0] + w[2]), 0.5 * (w[1] + w[3])
        if 215 < w[1] < 222 and w[0] < 300 and float(v).is_integer() \
                and 0 <= v <= 8:                       # panel (a) x ticks
            xs.append(cx)
            xv.append(v)
        elif w[2] <= 134.6 and 90 < w[1] < 215:        # panel (a) y ticks
            ys.append(cy)
            yv.append(v)
    px = np.polyfit(xs, xv, 1)
    py = np.polyfit(ys, yv, 1)
    rx = np.max(np.abs(np.polyval(px, xs) - xv))
    ry = np.max(np.abs(np.polyval(py, ys) - yv))
    print(f"calibration: {len(xs)} x-ticks residual {rx:.4f}, "
          f"{len(ys)} y-ticks residual {ry:.4f}")
    assert rx < TOL * 8 and ry < TOL, "tick calibration failed"
    return (lambda x: np.polyval(px, x)), (lambda y: np.polyval(py, y))


def main():
    p = fitz.open(PDF)[PAGE - 1]
    fx, fy = calib(p.get_text("words"))
    counts = {}
    for D in p.get_drawings():
        c = D.get("color")
        if c is None:
            continue
        key = tuple(round(v, 2) for v in c)
        if key not in COMP:
            continue
        pts = []
        for it in D["items"]:
            if it[0] == "l":
                pts.append((it[1].x, it[1].y))
                pts.append((it[2].x, it[2].y))
        pts = np.array(pts)
        e = fx(pts[:, 0])
        b = fy(pts[:, 1])
        m = (e >= -0.01) & (e <= EMAX)
        if m.sum() < 3:
            continue
        o = np.argsort(e[m])
        e, b = e[m][o], b[m][o]
        # collapse duplicate end-points
        keep = np.r_[True, np.diff(e) > 1e-4]
        e, b = e[keep], b[keep]
        name, note = COMP[key]
        k = counts.get(name, 0)
        counts[name] = k + 1
        fn = os.path.join(HERE, f"zp2013_{name}_path{k}.csv")
        with open(fn, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["e_total_strain", name, "their_component", note])
            for ee, bb in zip(e, b):
                w.writerow([f"{ee:.4f}", f"{bb:.4f}"])
        b05 = float(np.interp(0.5, e, b)) if e.max() >= 0.5 else np.nan
        print(f"{name} path{k}: {len(e):3d} pts, e in "
              f"[{e.min():.3f},{e.max():.3f}], b(e=0.5)={b05:+.3f}  "
              f"({note})")


if __name__ == "__main__":
    main()
