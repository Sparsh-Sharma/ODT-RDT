#!/usr/bin/env python3
"""Qualitative high-strain contrast of SC-ODT and clock-relaxed ODT, from
the real single-realization lines (hero_lines.npz, caro job 4461723).

The small-scale (high-pass) upwash structure is shown in the high-strain
window, in the stretched coordinate y/half-width so the compressed line
fills the frame. SC-ODT smooths to a few coarse bands; the clock-relaxed
line stays grainy, its sub-scale relaxation events keeping fine structure
alive. Qualitative -- amplitudes are normalised; the two are different
realisations, so this shows the clock's structural signature, not a
pointwise difference. NOT a main-paper figure.

  python fig_hero_diff_contour.py           colour
  python fig_hero_diff_contour.py --mono     greyscale
"""
import os
import sys

import numpy as np
from scipy.ndimage import uniform_filter1d

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from figstyle_jfm import FULL, plt, save  # noqa: E402

MONO = "--mono" in sys.argv
EMIN, EMAX, NX, WIN = 2.5, 4.0, 320, 25
D = np.load(os.path.join(HERE, "hero_lines.npz"))
yg = D["yg"]
CASES = [("hero_scodt", "SC-ODT", "smooths to a few coarse bands"),
         ("hero_clock", "clock-relaxed ODT",
          "stays grainy: sub-scale events\nkeep fine structure alive")]


def hp_field(key, comp="V"):
    t, yh, F = D["t_" + key], D["yh_" + key], D[comp + "_" + key]
    sel = np.where((t >= EMIN) & (t <= EMAX))[0]
    xi = np.linspace(-1, 1, NX)
    out = np.full((NX, len(sel)), np.nan)
    for k, j in enumerate(sel):
        col = F[:, j]
        m = np.isfinite(col)
        if m.sum() < 40:
            continue
        prof = np.interp(xi, yg[m] / yh[j], col[m], left=np.nan, right=np.nan)
        ins = np.isfinite(prof)
        out[ins, k] = prof[ins] - uniform_filter1d(prof[ins], WIN,
                                                    mode="nearest")
    return t[sel], xi, out


def main():
    cmap = "gray" if MONO else "RdBu_r"
    fig, axs = plt.subplots(1, 2, figsize=(FULL, 2.6), sharey=True)
    s = np.nanstd(np.concatenate(
        [hp_field(c[0])[2][np.isfinite(hp_field(c[0])[2])] for c in CASES]))
    for ax, (key, name, note) in zip(axs, CASES):
        te, xi, H = hp_field(key)
        ax.imshow(np.ma.masked_invalid(H), aspect="auto", origin="lower",
                  extent=[te[0], te[-1], -1, 1], cmap=cmap,
                  vmin=-2.2 * s, vmax=2.2 * s, interpolation="nearest",
                  rasterized=True)
        ax.set_title(name, fontsize=9, pad=4)
        ax.set_xlabel(r"total strain $e = St$")
        ax.set_xticks([2.5, 3.0, 3.5, 4.0])
        ax.text(0.5, -0.42, note, transform=ax.transAxes,
                ha="center", va="top", fontsize=6.8,
                color=("0.15" if MONO else "0.1"))
    axs[0].set_ylabel(r"$y\,/\,$half-width")
    axs[0].set_yticks([-1, 0, 1])
    fig.suptitle("fine-scale upwash structure at high strain (single line)",
                 fontsize=8.5, y=1.02)
    fig.tight_layout(pad=0.4, w_pad=1.4)
    fig.subplots_adjust(bottom=0.30, top=0.86)
    save(fig, os.path.join(HERE, "fig_hero_diff_contour"
                           + ("_mono" if MONO else "")))


if __name__ == "__main__":
    main()
