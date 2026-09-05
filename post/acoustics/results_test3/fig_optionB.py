#!/usr/bin/env python3
"""Option-B (unequal triplet-map images) summary figure: baseline vs
mapMidFrac 1/2 and 2/3.  Medians over 1024 realizations, 95% bootstrap CIs.
Writes fig_optionB.{png,pdf} next to this script."""
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = [("classic map (1/3,1/3,1/3)", "bands_homogeneousStrain2.npz", "k", "o"),
         ("images (1/4,1/2,1/4)", "bands_homogeneousStrain2B50.npz", "C0", "s"),
         ("images (1/6,2/3,1/6)", "bands_homogeneousStrain2B.npz", "C3", "^")]
NBOOT = 2000


def med_ci(x, rng):
    m = np.median(x)
    idx = rng.integers(0, x.size, (NBOOT, x.size))
    mb = np.median(x[idx], axis=1)
    return m, np.percentile(mb, 2.5), np.percentile(mb, 97.5)


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axs = plt.subplots(1, 4, figsize=(14.5, 3.7))
rng = np.random.default_rng(0)
for lab, fn, col, mk in CASES:
    d = np.load(os.path.join(HERE, fn))
    es = d["strains"]
    Alo = d["E2lo"] / d["Eplo"]
    Ahi = d["E2hi"] / d["Ephi"]
    series = [("$A_{low}$ ($k_2$ 30-100)", Alo),
              ("$A_{high}$ ($k_2$ 300-800)", Ahi),
              ("transmission $A_{high}/A_{low}$", Ahi / Alo),
              ("$u_2^2/2k_t$", d["u2frac"])]
    for ax, (title, arr) in zip(axs, series):
        m, lo, hi = np.array([med_ci(arr[:, j], rng) for j in range(len(es))]).T
        ax.errorbar(es, m, yerr=[m - lo, hi - m], fmt=mk + "-", color=col,
                    capsize=3, ms=4, lw=1.2, label=lab)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("total strain $e$")
axs[0].set_ylabel("median, 95% CI")
axs[2].axhline(1.0, color="k", lw=0.7, ls=":")
axs[3].axhline(1.0/3.0, color="k", lw=0.7, ls=":")
axs[3].axhline(0.6006, color="0.5", lw=0.7, ls="--")
axs[3].text(0.1, 0.605, "no-eddy LRR", fontsize=7, color="0.4")
axs[0].legend(fontsize=8, loc="upper left")
fig.suptitle("Option B (unequal triplet-map images), 1024 rlz/case", fontsize=11)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, f"fig_optionB.{ext}"), dpi=200)
print("saved fig_optionB.png/.pdf")
