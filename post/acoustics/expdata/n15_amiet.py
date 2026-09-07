"""Narayanan2015 Fig 4: digitise the Amiet prediction (red DASHED) curve.

Reuses the calibration from n15_extract.py. The solid measured curve is the
biggest red connected component; the dashes are the remaining red blobs.
Columns where a dash centre sits within a few px of the solid trace are
dropped as ambiguous (the two curves nearly touch at 2.5-4.5 kHz).
"""
import numpy as np, os
from PIL import Image
from scipy import ndimage

SCR = r"C:\Users\shar_sp\AppData\Local\Temp\claude\C--Users-shar-sp-Documents-GitHub-ODT-RDT\72078209-3184-44f4-936b-faab74533cac\scratchpad"
OUT = r"C:\Users\shar_sp\Documents\GitHub\ODT-RDT\post\acoustics\expdata"
a = np.array(Image.open(os.path.join(SCR, "n15_p13_img0.png")).convert("RGB")).astype(int)
R, G, B = a[..., 0], a[..., 1], a[..., 2]

XL, XR, YT, YB = 100.0, 694.0, 14.0, 374.0
def px2logf(x): return 2.0 + (x - XL) * 2.0 / (XR - XL)
def px2dB(y):  return 75.0 - (y - YT) * 50.0 / (YB - YT)

red = (R > 180) & (G < 90) & (B < 90)
lab, n = ndimage.label(red, structure=np.ones((3, 3), int))
sizes = ndimage.sum(red, lab, range(1, n + 1))
big = int(np.argmax(sizes)) + 1
solid = lab == big
dash = red & ~solid
print("red comps:", n, "solid px:", int(solid.sum()), "dash px:", int(dash.sum()))

# blob stats for the dashes
dlab, dn = ndimage.label(dash, structure=np.ones((3, 3), int))
dsz = ndimage.sum(dash, dlab, range(1, dn + 1))
print("dash blobs:", dn, "size min/med/max:", int(dsz.min()), int(np.median(dsz)), int(dsz.max()))
# drop tiny specks (<4 px) as anti-aliasing noise
keep = np.isin(dlab, np.where(dsz >= 4)[0] + 1)
dash = dash & keep

# solid trace per column (same tracker settings as n15_extract) for ambiguity test
def column_trace(mask, x0, x1, y_start, max_jump=8):
    ys, prev = {}, y_start
    for x in range(x0, x1 + 1):
        rows = np.where(mask[:, x])[0]
        if len(rows) == 0: continue
        runs, s = [], rows[0]
        for i in range(1, len(rows)):
            if rows[i] != rows[i - 1] + 1:
                runs.append((s, rows[i - 1])); s = rows[i]
        runs.append((s, rows[-1]))
        centers = [(r0 + r1) / 2 for r0, r1 in runs]
        c = min(centers, key=lambda v: abs(v - prev))
        if abs(c - prev) > max_jump and x > x0 + 3:
            continue
        ys[x] = c; prev = c
    return ys

ys_solid = column_trace(solid, 101, 693, 108, max_jump=8)

# dash per-column centres. XCUT: beyond x~590-600 (f ~4.5 kHz) the dashed curve merges
# into the thick solid measured band (verified in raw zoom) and the only remaining red
# fragments (x>=682) are debris of the spiky measured tangle -> ambiguous, dropped.
XCUT = 620
ys_dash, drop_near, drop_multi = {}, 0, 0
xs_all = []
for x in range(101, XCUT + 1):
    rows = np.where(dash[:, x])[0]
    if len(rows) == 0: continue
    runs, s = [], rows[0]
    for i in range(1, len(rows)):
        if rows[i] != rows[i - 1] + 1:
            runs.append((s, rows[i - 1])); s = rows[i]
    runs.append((s, rows[-1]))
    centers = [(r0 + r1) / 2 for r0, r1 in runs]
    if len(centers) > 1:
        # more than one dash run in a column: keep only if all within 6 px (thick dash), else drop
        if max(centers) - min(centers) > 6:
            drop_multi += 1; continue
        c = float(np.mean(centers))
    else:
        c = centers[0]
    # ambiguity: dash centre too close to solid trace -> can't be sure it's not bleed
    if x in ys_solid and abs(c - ys_solid[x]) < 4:
        drop_near += 1; continue
    ys_dash[x] = c
print("dash columns kept:", len(ys_dash), "dropped near-solid:", drop_near, "multi:", drop_multi)
xs = sorted(ys_dash)
print("dash x range:", xs[0], xs[-1], "-> f", 10**px2logf(xs[0]), 10**px2logf(xs[-1]))

def to_csv(ys, fname, nb, min_cols=2):
    xs = np.array(sorted(ys))
    logf = px2logf(xs)
    dB = px2dB(np.array([ys[x] for x in xs]))
    edges = np.linspace(logf.min(), logf.max(), nb + 1)
    fs, ds, sp = [], [], []
    for i in range(nb):
        m = (logf >= edges[i]) & (logf < edges[i + 1])
        if m.sum() < min_cols: continue
        fs.append(10 ** (0.5 * (edges[i] + edges[i + 1])))
        ds.append(np.median(dB[m]))
        sp.append(dB[m].max() - dB[m].min())
    with open(os.path.join(OUT, fname), "w") as f:
        f.write("f_Hz,SPL_dB\n")
        for xx, yy in zip(fs, ds):
            f.write(f"{xx:.1f},{yy:.2f}\n")
    print(f"{fname}: {len(fs)} pts, f {fs[0]:.0f}-{fs[-1]:.0f} Hz, med in-bin range {np.median(sp):.2f} dB, max {max(sp):.2f}")
    return np.array(fs), np.array(ds)

fA, dA = to_csv(ys_dash, "Narayanan2015_fig4_amiet60.csv", 55)

# offsets dashed - solid at 1, 2, 5 kHz (5 kHz is beyond the last dash: report anyway, flagged)
xs_s = np.array(sorted(ys_solid)); lf_s = px2logf(xs_s)
dB_s = px2dB(np.array([ys_solid[x] for x in xs_s]))
for f0 in (1000.0, 2000.0, 5000.0):
    lf0 = np.log10(f0)
    dsol = np.interp(lf0, lf_s, dB_s)
    dam = np.interp(lf0, np.log10(fA), dA)
    print(f"f={f0:.0f} Hz: Amiet {dam:.2f} dB, measured {dsol:.2f} dB, offset (Amiet-meas) {dam-dsol:+.2f} dB")

# check overlay: all three extracted curves on the raster
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7.67 * 1.6, 4.10 * 1.6), dpi=130)
ax.imshow(a.astype(np.uint8))
for fn, c in [("Narayanan2015_fig4_baseline60.csv", "blue"),
              ("Narayanan2015_fig4_background60.csv", "lime"),
              ("Narayanan2015_fig4_amiet60.csv", "black")]:
    d = np.loadtxt(os.path.join(OUT, fn), delimiter=",", skiprows=1)
    xsp = XL + (np.log10(d[:, 0]) - 2.0) / 2.0 * (XR - XL)
    yyp = YT + (75.0 - d[:, 1]) * (YB - YT) / 50.0
    ax.plot(xsp, yyp, "o", ms=3, color=c, mew=0)
ax.set_axis_off(); plt.tight_layout(pad=0)
plt.savefig(os.path.join(OUT, "Narayanan2015_fig4_amiet_check.png"))
print("check saved")
