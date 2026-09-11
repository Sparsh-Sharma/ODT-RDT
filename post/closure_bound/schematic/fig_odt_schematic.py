#!/usr/bin/env python3
"""Pedagogical schematic for sec:odt-review (request Sparsh 2026-09-16):
what an ODT eddy event is and how the model advances.

(a) an eddy event on a velocity profile: the triplet map on [y0, y0+l]
    (threefold compression, three images, central one reversed) plus the
    kernel addition c*K(y); exterior untouched.  Computed with the actual
    map, not drawn by hand.  Inset: the kernel K(y) = y - f(y).
(b) the two interleaved mechanisms in the space-time plane: instantaneous
    eddy events (Poisson in t, sampled sizes/positions) between which the
    fields advance by deterministic diffusion.
(c) the advancement cycle of the implementation (sample -> accept/reject
    -> map + kernels -> diffusion catch-up).

Monochrome per house rule 2026-09-08."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from figstyle_jfm import FULL, panel, plt, save  # noqa: E402

from matplotlib.patches import FancyArrowPatch, Rectangle  # noqa: E402


# ---------------------------------------------------------------- panel (a)
def triplet_f(yy, y0, l):
    """The triplet map pre-image point f(y) on [y0, y0+l]."""
    s = (yy - y0) / l
    out = np.where(s < 1.0 / 3.0, 3.0 * s,
                   np.where(s < 2.0 / 3.0, 2.0 - 3.0 * s, 3.0 * s - 2.0))
    return y0 + out * l


def make_profile(y, rng):
    u = np.zeros_like(y)
    for k, a in [(1, 0.55), (2, 0.42), (3, 0.30), (5, 0.16), (8, 0.07)]:
        u += a * np.sin(2 * np.pi * k * y + rng.uniform(0, 2 * np.pi))
    return u


def draw_panel_a(ax):
    rng = np.random.default_rng(7)
    y = np.linspace(0.0, 1.0, 4000)
    u = make_profile(y, rng)
    y0, l = 0.30, 0.42
    m = (y >= y0) & (y <= y0 + l)

    u_map = u.copy()
    u_map[m] = np.interp(triplet_f(y[m], y0, l), y, u)
    K = np.zeros_like(y)
    K[m] = y[m] - triplet_f(y[m], y0, l)
    c = 1.6                       # illustrative kernel amplitude
    off_b, off_a = 3.2, -0.9      # vertical offsets: before above, after below
    u_after = u_map + c * K

    ax.axvspan(y0, y0 + l, color="0.93", zorder=0)
    for x in (y0 + l / 3, y0 + 2 * l / 3):
        ax.axvline(x, color="0.78", lw=0.5, ls=(0, (2, 2)), zorder=1)
    ax.axvline(y0, color="0.55", lw=0.6, zorder=1)
    ax.axvline(y0 + l, color="0.55", lw=0.6, zorder=1)

    ax.plot(y, u + off_b, color="0.45", lw=1.0)
    ax.plot(y, u_after + off_a, color="k", lw=1.0)
    # the map alone (before the kernel addition), to separate the two steps
    ax.plot(y[m], u_map[m] + off_a, color="0.45", lw=0.6, ls=(0, (1, 1)))

    ax.annotate("", xy=(0.115, 0.10), xytext=(0.115, 1.55),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=0.8))
    ax.text(0.135, 0.85, "instantaneous\neddy event", fontsize=7,
            va="center")

    i05 = np.searchsorted(y, 0.06)
    ax.text(0.06, u[i05] + off_b + 0.75, "before", fontsize=7.5,
            color="0.35", ha="center")
    ax.text(0.06, u_after[i05] + off_a - 0.95, "after", fontsize=7.5,
            ha="center")
    ax.text(0.40, 4.95, r"$[y_0,\,y_0+l\,]$", fontsize=7,
            ha="right")
    for x, lab in [(y0 + l / 6, "1"), (y0 + l / 2, "2"),
                   (y0 + 5 * l / 6, "3")]:
        ax.text(x, -3.35, lab, fontsize=7, ha="center", color="0.3")

    ax.set_xlim(0, 1)
    ax.set_ylim(-3.9, 5.6)
    ax.set_xlabel(r"$y$")
    ax.set_ylabel(r"$u_i(y)$")
    ax.set_yticks([])
    ax.set_xticks([])

    # inset: the kernel K(y)
    axk = ax.inset_axes([0.68, 0.73, 0.27, 0.19])
    axk.plot(y[m], K[m], color="k", lw=0.8)
    axk.axhline(0, color="0.7", lw=0.4)
    axk.set_xticks([])
    axk.set_yticks([])
    for s in axk.spines.values():
        s.set_linewidth(0.5)
    axk.set_title(r"$K(y)=y-f(y)$", fontsize=7, pad=1.5)


# ---------------------------------------------------------------- panel (b)
def draw_panel_b(ax):
    # hand-placed events (t, y0, l): Poisson-looking, top strip kept clear
    events = [(0.07, 0.50, 0.28), (0.18, 0.18, 0.12), (0.28, 0.33, 0.22),
              (0.42, 0.20, 0.44), (0.545, 0.62, 0.18), (0.66, 0.25, 0.34),
              (0.78, 0.50, 0.10), (0.89, 0.18, 0.26)]
    for (te, y0, l) in events:
        ax.add_patch(Rectangle((te - 0.009, y0), 0.018, l,
                               facecolor="0.25", edgecolor="none", zorder=3))

    ax.text(0.03, 0.035, "between events:  "
            r"$\partial_t u_i=\nu\,\partial_y^2 u_i$",
            fontsize=7, ha="left", va="bottom", color="0.25")
    ax.annotate("instantaneous eddy events,\nPoisson in $t$, sampled "
                r"$(y_0,\,l\,)$",
                xy=(0.545, 0.805), xytext=(0.97, 0.985), fontsize=7,
                arrowprops=dict(arrowstyle="-", color="0.4", lw=0.5,
                                shrinkB=2),
                ha="right", va="top")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$y$")
    ax.set_xticks([])
    ax.set_yticks([])


# ---------------------------------------------------------------- panel (c)
def box(ax, x, y, w, h, text, fs=7, fc="white"):
    ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, facecolor=fc,
                           edgecolor="k", lw=0.6, zorder=2,
                           clip_on=False))
    ax.text(x, y, text, fontsize=fs, ha="center", va="center", zorder=3)


def draw_panel_c(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.0, 1.06, "(c)", fontstyle="italic", fontsize=8.5,
            va="top", ha="left")

    ys = [0.86, 0.60, 0.34, 0.08]
    box(ax, 0.5, ys[0], 0.92, 0.16,
        "sample next candidate time;\nadvance diffusion to it")
    box(ax, 0.5, ys[1], 0.92, 0.16,
        "draw candidate eddy $(y_0,\\,l\\,)$;\n"
        "timescale $\\tau(l\\,)$ from local field")
    box(ax, 0.5, ys[2], 0.92, 0.16,
        "accept with probability $\\propto 1/\\tau$\n"
        "(most candidates rejected)")
    box(ax, 0.5, ys[3], 0.92, 0.16,
        "apply event: triplet map $+\\,c_iK(y)$\n"
        "(momentum, energy conserving)", fc="0.92")

    for y1, y2 in zip(ys[:-1], ys[1:]):
        ax.add_patch(FancyArrowPatch((0.5, y1 - 0.08), (0.5, y2 + 0.08),
                                     arrowstyle="-|>", mutation_scale=7,
                                     color="k", lw=0.7, zorder=1))
    ax.text(0.545, (ys[2] + ys[3]) / 2, "accept", fontsize=7,
            ha="left", va="center", color="0.3")

    # reject: back to the top from the acceptance box, outside right
    ax.add_patch(FancyArrowPatch((0.965, ys[2]), (0.965, ys[0]),
                                 arrowstyle="-|>", mutation_scale=7,
                                 color="k", lw=0.7, zorder=1, clip_on=False,
                                 connectionstyle="arc3,rad=0.55"))
    ax.text(1.17, (ys[0] + ys[2]) / 2, "reject", fontsize=7, rotation=90,
            ha="center", va="center", color="0.3")
    # after an event: back to the top, outside left
    ax.add_patch(FancyArrowPatch((0.035, ys[3]), (0.035, ys[0]),
                                 arrowstyle="-|>", mutation_scale=7,
                                 color="k", lw=0.7, zorder=1, clip_on=False,
                                 connectionstyle="arc3,rad=-0.45"))


# ------------------------------------------------------------------- figure
def main():
    fig = plt.figure(figsize=(FULL, 2.55))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.20, 1.02],
                          wspace=0.34, left=0.045, right=0.955,
                          top=0.93, bottom=0.17)
    a = fig.add_subplot(gs[0])
    b = fig.add_subplot(gs[1])
    c = fig.add_subplot(gs[2])
    draw_panel_a(a)
    draw_panel_b(b)
    draw_panel_c(c)
    panel(a, "a")
    panel(b, "b")
    save(fig, os.path.join(HERE, "fig_odt_schematic"))


if __name__ == "__main__":
    main()
