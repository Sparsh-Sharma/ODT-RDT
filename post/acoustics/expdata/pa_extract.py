#!/usr/bin/env python3
"""Extract measured symbols from Paterson-Amiet CR-2733 fig 13 (p086.png,
300 dpi scan) by detecting the WHITE INTERIOR HOLES of the open plotting
symbols (squares/triangles/diamonds/circles) - robust to symbols touching
the theory curves and to the little background flags.

Classification of the hole shape:
  circle   : 4-fold symmetric, filled corner fraction high, round
  square   : corners filled (extent ~1 in bbox), aspect ~1
  diamond  : rotated square: extent ~0.5, aspect ~1
  tri_up   : centroid BELOW bbox centre (mass at bottom), extent ~0.5
  tri_dn   : centroid ABOVE bbox centre
Outputs PatersonAmiet_fig13_<label>.csv (f_Hz, SPL_dB) per series +
overlay proof pa76_png/fig13_check.png.
"""
import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage

X0, PXHZ = 406.0, (2000.0 - 406.0) / 5000.0
Y0, PXDB = 505.0, (2232.0 - 505.0) / 70.0
im = Image.open("pa76_png/p086.png").convert("L")
a = np.asarray(im)
# plot interior only
box = (slice(int(Y0) + 2, 2232 - 2), slice(int(X0) + 2, 2000 - 2))
dark = a[box] < 128
# white regions NOT connected to the plot border = candidate holes
white = ~dark
lab, n = ndimage.label(white)
border = np.unique(np.r_[lab[0, :], lab[-1, :], lab[:, 0], lab[:, -1]])
objs = ndimage.find_objects(lab)
rows = []
for i in range(1, n + 1):
    if i in border:
        continue
    sl = objs[i - 1]
    h = sl[0].stop - sl[0].start
    w = sl[1].stop - sl[1].start
    npx = int((lab[sl] == i).sum())
    if not (8 <= h <= 45 and 8 <= w <= 45 and npx >= 30):
        continue
    m = lab[sl] == i
    ys, xs = np.nonzero(m)
    cy, cx = ys.mean(), xs.mean()
    extent = npx / (h * w)
    # vertical mass asymmetry: (centroid - centre)/h
    va = (cy - (h - 1) / 2.0) / h
    # corner fill: fraction of the 4 bbox corner quarters occupied
    q = h // 3, w // 3
    corners = (m[: q[0], : q[1]].mean() + m[: q[0], -q[1]:].mean()
               + m[-q[0]:, : q[1]].mean() + m[-q[0]:, -q[1]:].mean()) / 4
    if extent > 0.72 and corners > 0.45:
        kind = "square"
    elif va > 0.06:
        kind = "tri_up"
    elif va < -0.06:
        kind = "tri_dn"
    elif extent > 0.62:
        kind = "circle"
    else:
        kind = "diamond"
    # page coords
    py = box[0].start + sl[0].start + cy
    px = box[1].start + sl[1].start + cx
    f = (px - X0) / PXHZ
    db = 90.0 - (py - Y0) / PXDB
    rows.append((kind, f, db, px, py, extent, va, corners, h, w))

# drop detections inside the legend and the in-plot text labels
EXCL = [(1585, 480, 1835, 1060),    # legend block + "VELOCITY, m/sec"
        (920, 520, 1540, 675),      # "phiM = 90 DEG" / "THEORY"
        (1390, 1200, 1485, 1270),   # "120" label
        (1190, 1265, 1275, 1340),   # "90" label
        (1165, 1580, 1250, 1655),   # "60" label
        (985, 1780, 1070, 1855),    # "40" label
        (1715, 1105, 1890, 1175)]   # "165 m/sec" label
rows = [r for r in rows
        if not any(x0 <= r[3] <= x1 and y0 <= r[4] <= y1
                   for x0, y0, x1, y1 in EXCL)]

# The five series occupy vertically disjoint bands at every frequency, so
# classify purely by position between hand-read separator polylines (the
# 1976 scan is too noisy for reliable hole-shape moments).
SEPS = [  # from top: above s0 -> U165, between s0,s1 -> U120, ...
    ([150, 600, 1000, 1600, 2000, 2600, 3200, 4500],
     [77, 78.5, 76, 70.5, 67, 65, 62, 58]),
    ([150, 500, 1000, 1500, 2000, 2600, 3200, 3800],
     [72, 73.5, 68, 62, 58.5, 55, 53.5, 52]),
    ([150, 400, 800, 1200, 1600, 2000, 2600, 3800],
     [66, 67, 62.5, 54, 49, 46.5, 44, 40]),
    ([150, 400, 800, 1200, 1500, 2100],
     [56, 57.5, 50.75, 41, 36, 32]),
]
ORDER = ["U165", "U120", "U90", "U60", "U40"]


def classify(f, db):
    for i, (fs, dbs) in enumerate(SEPS):
        if db > np.interp(f, fs, dbs):
            return ORDER[i]
    return ORDER[-1]


rows = [(classify(f, db), f, db, px, py, *rest)
        for _k, f, db, px, py, *rest in rows]

# merge duplicate detections of one symbol (broken outlines -> two holes)
merged = []
for r in sorted(rows, key=lambda r: (r[0], r[1])):
    if merged and merged[-1][0] == r[0] \
            and abs(merged[-1][3] - r[3]) < 14 \
            and abs(merged[-1][4] - r[4]) < 14:
        continue
    merged.append(r)
rows = merged

series = {"U165": "U165", "U120": "U120", "U90": "U90",
          "U60": "U60", "U40": "U40"}
proof = Image.open("pa76_png/p086.png").convert("RGB")
dr = ImageDraw.Draw(proof)
cols = {"U165": (230, 0, 0), "U120": (0, 150, 0), "U90": (0, 0, 230),
        "U60": (200, 120, 0), "U40": (170, 0, 170)}
for kind, lbl in series.items():
    pts = sorted([(f, db) for k, f, db, *_ in rows if k == kind])
    with open(f"PatersonAmiet_fig13_{lbl}.csv", "w") as fh:
        fh.write("f_Hz,SPL_dB\n")
        for f, db in pts:
            fh.write(f"{f:.0f},{db:.2f}\n")
    print(lbl, len(pts), "points",
          f"f {pts[0][0]:.0f}-{pts[-1][0]:.0f} Hz" if pts else "")
for kind, f, db, px, py, *_ in rows:
    dr.ellipse([px - 6, py - 6, px + 6, py + 6], outline=cols[kind], width=2)
proof.save("pa76_png/fig13_check.png")
print("proof: pa76_png/fig13_check.png; total symbols:", len(rows))
