#!/usr/bin/env python3
"""Space-time propagation of the ODT line as a flowing ridgeline of the
actual velocity profile, one column per model variant (real runs, same
seed 22, eddies on; caro job 4461723 -> hero_lines.npz).

Each stacked curve is the line's velocity profile at one instant; time
runs upward; colour flows with time. Standard ODT keeps its full width;
SC-ODT and clock-relaxed taper as the line physically compresses under
strain. Per-profile normalisation keeps the eddy structure visible as the
fluctuations decay.

  python fig_hero_ribbon.py           colour (proposals)
  python fig_hero_ribbon.py --mono    greyscale (paper)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from figstyle_jfm import FULL, plt, save  # noqa: E402

MONO = "--mono" in sys.argv
COMP = "V"
NPROF = 55
STEP = 1.0
AMP = 1.35
D = np.load(os.path.join(HERE, "hero_lines.npz"))
yg = D["yg"]
CASES = [("hero_std", "standard ODT"),
         ("hero_scodt", "SC-ODT"),
         ("hero_clock", "clock-relaxed ODT")]
CMAP = plt.get_cmap("Greys") if MONO else plt.get_cmap("viridis")


def main():
    fig, axs = plt.subplots(1, 3, figsize=(FULL, 3.5), sharey=True)
    for ax, (key, name) in zip(axs, CASES):
        F = D[COMP + "_" + key]
        yh = D["yh_" + key]
        nt = F.shape[1]
        sel = np.linspace(0, nt - 1, NPROF).astype(int)
        for k, j in enumerate(sel):
            prof = F[:, j]
            m = np.isfinite(prof)
            if m.sum() < 5:
                continue
            y = yg[m]
            p = prof[m]
            p = p / (np.nanmax(np.abs(p)) + 1e-9)
            base = k * STEP
            curve = base + AMP * p
            col = CMAP(0.12 + 0.8 * k / (NPROF - 1))
            ax.fill_between(y, base, curve, color=col, alpha=0.9,
                            lw=0, zorder=k)
            ax.plot(y, curve, color="k", lw=0.25, alpha=0.6, zorder=k)
        # compression cone (physical half-width vs time index)
        kk = np.linspace(0, NPROF - 1, NPROF)
        yhs = yh[sel]
        ax.plot(yhs, kk * STEP, color="k", lw=0.8, zorder=NPROF + 1)
        ax.plot(-yhs, kk * STEP, color="k", lw=0.8, zorder=NPROF + 1)
        ax.set_xlim(-0.55, 0.55)
        ax.set_title(name, fontsize=9, pad=5)
        ax.set_xlabel(r"line coordinate $y$")
        ax.set_xticks([-0.4, 0, 0.4])
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.set_yticks([])
    # strain axis on the first panel (e = St runs bottom->top)
    etick = [0, 1, 2, 3, 4]
    kpos = [E / 4.0 * (NPROF - 1) * STEP for E in etick]
    axs[0].spines["left"].set_visible(True)
    axs[0].spines["left"].set_bounds(0, (NPROF - 1) * STEP)
    axs[0].set_yticks(kpos)
    axs[0].set_yticklabels([str(E) for E in etick])
    axs[0].set_ylabel(r"total strain $e = St$")
    fig.tight_layout(pad=0.5, w_pad=1.0)
    save(fig, os.path.join(HERE, "fig_hero_ribbon"
                           + ("_mono" if MONO else "")))


if __name__ == "__main__":
    main()
