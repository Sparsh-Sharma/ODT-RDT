#!/usr/bin/env python3
"""Space-time propagation of the ODT line for three model variants, from
real single-realization runs (same seed 22, same initial field, eddies on;
caro job 4461723, reduced to hero_lines.npz):

  standard ODT        -- strain off: homogeneous turbulence, full-width line
  SC-ODT              -- plane strain + dilatation: the line compresses ~7x
  clock-relaxed ODT   -- SC-ODT plus the independent relaxation clock

Each strip is the upwash velocity along the line against time; the field is
normalised per time so the eddy structure stays visible as the fluctuations
decay, and the black envelope is the physical line half-width (the
compression under strain).

  python fig_hero_lines.py           monochrome (paper)
  python fig_hero_lines.py --color   colour (proposals)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from figstyle_jfm import FULL, plt, save  # noqa: E402
import matplotlib.colors as mcolors       # noqa: E402

COLOR = "--color" in sys.argv
COMP = "V"          # upwash u2 (V); "U" streamwise, "W" spanwise
D = np.load(os.path.join(HERE, "hero_lines.npz"))
yg = D["yg"]

CASES = [("hero_std", "standard ODT", "strain off"),
         ("hero_scodt", "SC-ODT", "plane strain, dilatation"),
         ("hero_clock", "clock-relaxed ODT", "strain + relaxation clock")]

CMAP = "RdBu" if COLOR else "gray"   # colour: blue positive, red negative


def colnorm(A):
    s = np.nanstd(A, axis=0, keepdims=True)
    s[s == 0] = np.nan
    return A / s


def main():
    n = len(CASES)
    h = 0.92 * n + 0.6
    fig, axs = plt.subplots(n, 1, figsize=(FULL, h), sharex=True)
    cmap = plt.get_cmap(CMAP).copy()
    cmap.set_bad(color=(1, 1, 1, 0))          # outside the line: transparent
    vlim = 2.6
    for ax, (key, name, sub) in zip(axs, CASES):
        t = D["t_" + key]
        yh = D["yh_" + key]
        F = colnorm(D[COMP + "_" + key])       # (Ny, Nt)
        te = np.append(t, 2 * t[-1] - t[-2])   # cell edges in time
        ye = np.append(yg, 2 * yg[-1] - yg[-2])
        pcm = ax.pcolormesh(te, ye, np.ma.masked_invalid(F), cmap=cmap,
                            vmin=-vlim, vmax=vlim, rasterized=True,
                            shading="flat")
        # physical line envelope (the compression under strain)
        ax.plot(t, yh, color="k", lw=0.7)
        ax.plot(t, -yh, color="k", lw=0.7)
        ax.set_ylim(-0.52, 0.52)
        ax.set_xlim(t[0], t[-1])
        ax.set_yticks([-0.4, 0, 0.4])
        ax.set_ylabel(r"$y$")
        ax.text(0.012, 0.9, name, transform=ax.transAxes, va="top",
                ha="left", fontsize=8.5, fontweight="bold",
                bbox=dict(boxstyle="square,pad=0.18", fc="white",
                          ec="0.6", lw=0.4, alpha=0.9))
        ax.text(0.012, 0.10, sub, transform=ax.transAxes, va="bottom",
                ha="left", fontsize=6.8, color="0.25",
                bbox=dict(boxstyle="square,pad=0.12", fc="white",
                          ec="none", alpha=0.75))
    # compression annotation on the SC-ODT strip
    axs[1].annotate("line compresses\nunder strain",
                    xy=(3.9, D["yh_hero_scodt"][-1]), xytext=(2.6, 0.34),
                    fontsize=6.6, ha="center", color="0.15",
                    arrowprops=dict(arrowstyle="-|>", color="0.3", lw=0.7))
    axs[-1].set_xlabel(r"total strain $e=St$   (time $t$ for standard ODT)")
    fig.tight_layout(pad=0.4, h_pad=0.5)
    cb = fig.colorbar(pcm, ax=list(axs), location="right", fraction=0.022,
                      pad=0.015, aspect=42, ticks=[-vlim, 0, vlim])
    cb.set_label(r"upwash velocity $u_2$ (normalised per time)", fontsize=7)
    cb.ax.set_yticklabels([r"$-$", "0", r"$+$"])
    cb.ax.tick_params(labelsize=7, length=2)
    save(fig, os.path.join(HERE, "fig_hero_lines_color" if COLOR
                           else "fig_hero_lines"))


if __name__ == "__main__":
    main()
