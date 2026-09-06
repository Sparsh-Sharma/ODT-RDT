#!/usr/bin/env python3
"""Hierarchical-kernel (Kerstein 2026-09-06) summary: baseline vs
nSubKernelLevels = 1, 2, 3.  Medians over 1024 realizations, 95% bootstrap
CIs; all four ensembles same (post-dilatation-fix) binary, paired seeds.
Writes fig_klevels.{png,pdf} and klevels_table.txt."""
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NB = os.path.join(HERE, "newbands")
CASES = [("standard ODT", "bands_homogeneousStrain2.npz", "k", "o"),
         ("level 1 (3 images)", "bands_homogeneousStrain2K1.npz", "C0", "s"),
         ("level 2 (+9)", "bands_homogeneousStrain2K2.npz", "C2", "^"),
         ("level 3 (+27)", "bands_homogeneousStrain2K3.npz", "C3", "d")]
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
lines = ["Hierarchical kernel levels, medians [95% CI], 1024 rlz each",
         f"{'case':>20} {'A_low(1)':>16} {'A_high(1)':>16} {'A_high(3)':>16} "
         f"{'T(3.9)':>16} {'u2/2kt(3.9)':>12}"]
for lab, fn, col, mk in CASES:
    d = np.load(os.path.join(NB, fn))
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
    j1 = int(np.argmin(np.abs(es - 1.0)))
    j3 = int(np.argmin(np.abs(es - 3.0)))
    j4 = len(es) - 1
    f = lambda x: f"{np.median(x):.2f}"
    lines.append(f"{lab:>20} {f(Alo[:, j1]):>16} {f(Ahi[:, j1]):>16} "
                 f"{f(Ahi[:, j3]):>16} {f(Ahi[:, j4]/Alo[:, j4]):>16} "
                 f"{np.median(d['u2frac'][:, j4]):>12.3f}")
axs[0].set_ylabel("median, 95% CI")
axs[2].axhline(1.0, color="k", lw=0.7, ls=":")
axs[0].legend(fontsize=8, loc="upper left")
fig.suptitle("Hierarchical kernel isotropization (levels 0-3), 1024 rlz/case",
             fontsize=11)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, f"fig_klevels.{ext}"), dpi=200)
txt = "\n".join(lines)
print(txt)
open(os.path.join(HERE, "klevels_table.txt"), "w").write(txt + "\n")
print("saved fig_klevels.png/.pdf")
