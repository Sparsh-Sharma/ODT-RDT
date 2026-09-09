#!/usr/bin/env python3
"""Digitise Zusi & Perot 2014 (Phys. Fluids 26, 115103), FIG. 7, p. 10:
normalised anisotropy tensor against time for axisymmetric contraction
(panel a, AXC IC1) and expansion (panel b, AXE IC2) at low/medium/high
strain rates.  Vector figure; one PDF path per colour (= strain rate)
per panel, each path holding all three components:
  b11 solid  -> one very long connected chain,
  b22 dashed -> chains of ~100-270 segments (each dash ~10 pt),
  b33 dotted -> chains of 10-99 segments (each dot ~1 pt).
Straining runs from t = 12.0 until 12.8 (high), 15.2 (medium), 28.0
(low) s; the colour -> rate assignment is identified from the time of
the b11 extremum and printed, then checked against those end times.

Their strain measure: e_ZP = S (t - 12) with S = Sk0/eps0 / T0 from their
Table I (S = the maximum absolute mean strain rate = the axial rate).
Every case ends at e_ZP = 0.5.  The CSVs carry t, e_ZP and b in THEIR
units; any conversion to the manuscript's strain measure is applied at
overlay time, not here.

Usage: python extract_zp2014.py   (writes zp2014_<AXC|AXE>_<rate>_<bii>.csv)
"""
import csv
import os

import fitz
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "ZusiPerot2014.pdf")
PAGE = 10
TOL = 0.01
COLS = [(0.93, 0.12, 0.14), (0.22, 0.32, 0.64), (0.18, 0.6, 0.4)]
T_END = {"high": 12.8, "med": 15.2, "low": 28.0}
# Table I: (T0, {rate: Sk0/eps0})
PANEL = {"a": ("AXC", 5.388, {"low": 0.168, "med": 0.842, "high": 3.367}),
         "b": ("AXE", 5.252, {"low": 0.164, "med": 0.821, "high": 3.282})}


def calib(words, panel):
    xs, xv, ys, yv = [], [], [], []
    for w in words:
        s = w[4].replace("−", "-").replace("–", "-")
        try:
            v = float(s)
        except ValueError:
            continue
        cx, cy = 0.5 * (w[0] + w[2]), 0.5 * (w[1] + w[3])
        inpan = (w[0] < 305) if panel == "a" else (w[0] > 325)
        if 219 < w[1] < 222 and inpan and 10 <= v <= 100:
            xs.append(cx)
            xv.append(v)
        yedge = 131.8 if panel == "a" else 331.6
        if abs(w[2] - yedge) < 0.5 and 90 < w[1] < 220 and -0.2 <= v <= 0.15:
            ys.append(cy)
            yv.append(v)
    px, py = np.polyfit(xs, xv, 1), np.polyfit(ys, yv, 1)
    rx = np.max(np.abs(np.polyval(px, xs) - xv))
    ry = np.max(np.abs(np.polyval(py, ys) - yv))
    print(f"panel {panel}: {len(xs)} x-ticks residual {rx:.3f} s, "
          f"{len(ys)} y-ticks residual {ry:.4f}")
    assert rx < 0.5 and ry < TOL, "tick calibration failed"
    return (lambda x: np.polyval(px, x)), (lambda y: np.polyval(py, y))


def chains(D):
    S = np.array([(it[1].x, it[1].y, it[2].x, it[2].y)
                  for it in D["items"] if it[0] == "l"])
    conn = np.hypot(S[1:, 0] - S[:-1, 2], S[1:, 1] - S[:-1, 3]) < 1e-3
    br = np.where(~conn)[0]
    starts, ends = np.r_[0, br + 1], np.r_[br, len(S) - 1]
    return [S[a:b + 1] for a, b in zip(starts, ends)]


def classify(ch):
    n = np.array([len(c) for c in ch])
    solid = [ch[int(np.argmax(n))]]
    dashed = [c for c, k in zip(ch, n) if k >= 100 and c is not solid[0]]
    dotted = [c for c, k in zip(ch, n) if 10 <= k < 100]
    return {"b11": solid, "b22": dashed, "b33": dotted}


def to_xy(segs, fx, fy):
    P = np.concatenate([np.r_[c[:, :2], c[:, 2:]] for c in segs])
    t, b = fx(P[:, 0]), fy(P[:, 1])
    o = np.argsort(t)
    t, b = t[o], b[o]
    keep = np.r_[True, np.diff(t) > 1e-3]
    return t[keep], b[keep]


def main():
    p = fitz.open(PDF)[PAGE - 1]
    words = p.get_text("words")
    for panel, (name, T0, rates) in PANEL.items():
        fx, fy = calib(words, panel)
        for D in p.get_drawings():
            c = D.get("color")
            if c is None or tuple(round(v, 2) for v in c) not in COLS:
                continue
            if (D["rect"].x0 < 300) != (panel == "a"):
                continue
            key = tuple(round(v, 2) for v in c)
            comp = classify(chains(D))
            curves = {k: to_xy(v, fx, fy) for k, v in comp.items() if v}
            # rate from the b11 extremum time (AXC: minimum; AXE: maximum)
            t11, b11 = curves["b11"]
            m = (t11 >= 12) & (t11 <= 40)
            j = np.argmin(b11[m]) if name == "AXC" else np.argmax(b11[m])
            text = t11[m][j]
            rate = min(T_END, key=lambda r: abs(T_END[r] - text))
            S = rates[rate] / T0
            print(f"{name} colour {key}: b11 extremum at t={text:.2f} -> "
                  f"'{rate}' (t_end {T_END[rate]}), S={S:.4f} s^-1, "
                  f"e_ZP(t_end)={S * (T_END[rate] - 12):.3f}")
            for comp_name, (t, b) in curves.items():
                w = (t >= 12.0) & (t <= T_END[rate] + 0.02)
                fn = os.path.join(HERE, f"zp2014_{name}_{rate}_{comp_name}.csv")
                with open(fn, "w", newline="") as f:
                    wr = csv.writer(f)
                    wr.writerow(["t_s", "e_zp", comp_name,
                                 f"{name} IC{'1' if name == 'AXC' else '2'} "
                                 f"Sk0/eps0={rates[rate]} S={S:.4f}"])
                    for tt, bb in zip(t[w], b[w]):
                        wr.writerow([f"{tt:.4f}", f"{S * (tt - 12):.4f}",
                                     f"{bb:.4f}"])
                print(f"    {comp_name}: {w.sum():4d} pts in strain window, "
                      f"b(t_end)={np.interp(T_END[rate], t, b):+.3f}")


if __name__ == "__main__":
    main()
