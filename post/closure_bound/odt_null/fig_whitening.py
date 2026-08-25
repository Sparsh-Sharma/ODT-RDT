"""Figure for the null-test note: b_ii(t) for the three IC-whitening variants.

Three panels sharing one y-axis; fixed component colors (validated palette),
distinct markers as secondary encoding, direct labels, jackknife error bars.
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(HERE, "whitening_bt.npz"))
times = d["times"]

COL = ["#2a78d6", "#eb6834", "#1baf7a"]        # u, v, w (fixed order)
MRK = ["o", "s", "^"]
LBL = [r"$b_{11}$ ($u$)", r"$b_{22}$ ($v$, line)", r"$b_{33}$ ($w$)"]
PANELS = [("forward", "ordered Cholesky (u,v,w)\n(original code)"),
          ("reversed", "ordered Cholesky (w,v,u)\n(attribution test)"),
          ("symmetric", r"symmetric $C^{-1/2}$" + "\n(fixed code)")]

fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), sharey=True)
for ax, (key, title) in zip(axes, PANELS):
    b, se = d[key + "_b"], d[key + "_se"]
    ax.axhline(0.0, color="0.75", lw=0.8, zorder=0)
    for c in range(3):
        ax.errorbar(times[1:], b[1:, c], yerr=se[1:, c], color=COL[c],
                    marker=MRK[c], ms=4.5, lw=1.6, capsize=2.5,
                    label=LBL[c])
    ax.set_title(title, fontsize=9.5)
    ax.set_xlabel(r"$t$ (eddy time units)", fontsize=9)
    ax.tick_params(labelsize=8.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.92", lw=0.6, zorder=0)
axes[0].set_ylabel(r"component anisotropy $b_{ii}$", fontsize=9.5)
axes[0].annotate("last-whitened component\nover-energized",
                 xy=(0.62, 0.0146), xytext=(1.55, 0.0165), fontsize=8,
                 color="0.25", va="center",
                 arrowprops=dict(arrowstyle="-", color="0.55", lw=0.8))
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, ncol=3, fontsize=9, frameon=False,
           loc="upper center", bbox_to_anchor=(0.5, 1.02))
fig.suptitle("ODT nominal-HIT null test: the $b_{ii}$ artifact follows the "
             "IC whitening order (1024 realizations each)", fontsize=10.5,
             y=1.10)
fig.tight_layout(rect=(0, 0, 1, 0.97))
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, "fig_whitening_bt." + ext),
                bbox_inches="tight", dpi=180)
print("saved fig_whitening_bt.pdf/.png")
