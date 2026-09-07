"""Narayanan2015 Fig 4: digitise measured baseline (red solid) + background noise (orange dashed)."""
import numpy as np, os
from PIL import Image
from scipy import ndimage

SCR = r"C:\Users\shar_sp\AppData\Local\Temp\claude\C--Users-shar-sp-Documents-GitHub-ODT-RDT\72078209-3184-44f4-936b-faab74533cac\scratchpad"
OUT = r"C:\Users\shar_sp\Documents\GitHub\ODT-RDT\post\acoustics\expdata"
a = np.array(Image.open(os.path.join(SCR, "n15_p13_img0.png")).convert("RGB")).astype(int)
R, G, B = a[..., 0], a[..., 1], a[..., 2]

# calibration: x=100 -> 10^2 Hz, x=694 -> 10^4 Hz; y=14 -> 75 dB, y=374 -> 25 dB
XL, XR, YT, YB = 100.0, 694.0, 14.0, 374.0
def px2logf(x): return 2.0 + (x - XL) * 2.0 / (XR - XL)
def px2dB(y):  return 75.0 - (y - YT) * 50.0 / (YB - YT)

red = (R > 180) & (G < 90) & (B < 90)
orange = (R > 170) & (R < 255) & (G > 90) & (G < 170) & (B < 70)

# --- red: separate solid (one big connected component) from dashed ---
lab, n = ndimage.label(red, structure=np.ones((3, 3), int))
sizes = ndimage.sum(red, lab, range(1, n + 1))
big = int(np.argmax(sizes)) + 1
solid = lab == big
print("red components:", n, "biggest size", int(sizes.max()), "total red", int(red.sum()))

def column_trace(mask, x0, x1, y_start, max_jump=8):
    ys = {}
    prev = y_start
    for x in range(x0, x1 + 1):
        rows = np.where(mask[:, x])[0]
        if len(rows) == 0:
            continue
        # runs
        runs, s = [], rows[0]
        for i in range(1, len(rows)):
            if rows[i] != rows[i - 1] + 1:
                runs.append((s, rows[i - 1])); s = rows[i]
        runs.append((s, rows[-1]))
        centers = [(r0 + r1) / 2 for r0, r1 in runs]
        c = min(centers, key=lambda v: abs(v - prev))
        if abs(c - prev) > max_jump and x > x0 + 3:
            continue
        ys[x] = c
        prev = c
    return ys

# solid baseline starts at left edge ~62 dB -> y ~ 107
ys_solid = column_trace(solid, 101, 693, 108, max_jump=8)
print("solid columns:", len(ys_solid))

# --- orange dashed: per-column centers (single curve) ---
ys_bg = {}
for x in range(101, 694):
    rows = np.where(orange[:, x])[0]
    if len(rows) == 0: continue
    if rows.max() - rows.min() > 12:  # ambiguous (shouldn't happen)
        continue
    ys_bg[x] = rows.mean()
print("background columns:", len(ys_bg))

def to_csv(ys, fname, nb):
    xs = np.array(sorted(ys))
    logf = px2logf(xs)
    dB = px2dB(np.array([ys[x] for x in xs]))
    edges = np.linspace(logf.min(), logf.max(), nb + 1)
    fs, ds, sp = [], [], []
    for i in range(nb):
        m = (logf >= edges[i]) & (logf < edges[i + 1])
        if m.sum() < 2: continue
        fs.append(10 ** (0.5 * (edges[i] + edges[i + 1])))
        ds.append(np.median(dB[m]))
        sp.append(dB[m].max() - dB[m].min())
    with open(os.path.join(OUT, fname), "w") as f:
        f.write("f_Hz,SPL_dB\n")
        for xx, yy in zip(fs, ds):
            f.write(f"{xx:.1f},{yy:.2f}\n")
    print(f"{fname}: {len(fs)} pts, f {fs[0]:.0f}-{fs[-1]:.0f} Hz, med in-bin range {np.median(sp):.2f} dB")

to_csv(ys_solid, "Narayanan2015_fig4_baseline60.csv", 70)
to_csv(ys_bg, "Narayanan2015_fig4_background60.csv", 55)

# check plot
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7.67 * 1.6, 4.10 * 1.6), dpi=130)
ax.imshow(a.astype(np.uint8))
for fn, c in [("Narayanan2015_fig4_baseline60.csv", "blue"),
              ("Narayanan2015_fig4_background60.csv", "lime")]:
    d = np.loadtxt(os.path.join(OUT, fn), delimiter=",", skiprows=1)
    xs = XL + (np.log10(d[:, 0]) - 2.0) / 2.0 * (XR - XL)
    yy = YT + (75.0 - d[:, 1]) * (YB - YT) / 50.0
    ax.plot(xs, yy, "o", ms=3, color=c, mew=0)
ax.set_axis_off(); plt.tight_layout(pad=0)
plt.savefig(os.path.join(OUT, "Narayanan2015_fig4_check.png"))
print("check saved")
