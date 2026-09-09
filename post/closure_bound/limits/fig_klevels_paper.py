#!/usr/bin/env python3
"""Manuscript figure for sec:limits (main_v3.tex): hierarchical sub-scale
kernel relaxation — standard ODT vs levels 1-3, adaptive depth, and the
iterated variant.  Medians over 1024 paired-seed realizations, 95%
bootstrap CIs, one binary.  Data: newbands/*.npz (mirrored from
LEN_Extension post/acoustics/results_test3/newbands, campaigns of
2026-09-06).  Output: fig_klevels.{pdf,png} at true \\textwidth."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from figstyle_jfm import FULL, panel, plt, save  # noqa: E402

NB = os.path.join(HERE, "newbands")
# monochrome: baseline black solid; the fixed-depth 1/2/3-level family a
# light->dark grey ramp (deeper = darker); the other four each a distinct
# black line style.  (label, file, marker, colour, linestyle)
CASES = [
    ("SC-ODT (baseline)", "bands_homogeneousStrain2.npz", "o", "k", "-"),
    ("1 level (thirds)", "bands_homogeneousStrain2K1.npz", "s", "0.62",
     "-"),
    ("2 levels", "bands_homogeneousStrain2K2.npz", "^", "0.42", "-"),
    ("3 levels", "bands_homogeneousStrain2K3.npz", "D", "0.18", "-"),
    ("depth tied to eddy size", "bands_homogeneousStrain2KS.npz",
     "v", "k", (0, (6, 2))),
    ("2 levels, iterated 3x", "bands_homogeneousStrain2K2I3.npz",
     "P", "k", (0, (1, 1.4))),
    ("concurrent, $N{=}3$", "bands_homogeneousStrain2CR3.npz",
     "X", "k", (0, (5, 1.6, 1, 1.6))),
    ("concurrent + depth", "bands_homogeneousStrain2CRS.npz", "*", "k",
     (0, (4, 1.4, 1, 1.4, 1, 1.4))),
]
NBOOT = 2000


def med_ci(x, rng):
    m = np.median(x)
    idx = rng.integers(0, x.size, (NBOOT, x.size))
    mb = np.median(x[idx], axis=1)
    return m, np.percentile(mb, 2.5), np.percentile(mb, 97.5)


def main():
    rng = np.random.default_rng(0)
    fig, axs = plt.subplots(2, 2, figsize=(FULL, 3.9), sharex=True)
    axs = axs.ravel()
    for lab, fn, mk, col, ls in CASES:
        d = np.load(os.path.join(NB, fn))
        es = d["strains"]
        Alo = d["E2lo"] / d["Eplo"]
        Ahi = d["E2hi"] / d["Ephi"]
        series = [(r"$\mathcal{A}_{\rm lo}$", Alo),
                  (r"$\mathcal{A}_{\rm hi}$", Ahi),
                  (r"$\mathcal{A}_{\rm hi}/\mathcal{A}_{\rm lo}$", Ahi / Alo),
                  (r"$u_2^2/2k_t$", d["u2frac"])]
        for ax, (ylab, arr) in zip(axs, series):
            m, lo, hi = np.array([med_ci(arr[:, j], rng)
                                  for j in range(len(es))]).T
            ax.errorbar(es, m, yerr=[m - lo, hi - m], fmt=mk, ls=ls,
                        color=col, mfc=col, capsize=1.4, capthick=0.5,
                        elinewidth=0.5, ms=2.8, lw=0.9, label=lab)
            ax.set_ylabel(ylab)
    axs[2].axhline(1.0, color="k", lw=0.6, ls=":")
    for j, ax in enumerate(axs):
        if j >= 2:
            ax.set_xlabel("total strain $e$")
        panel(ax, "abcd"[j])
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, ncol=4, loc="lower center",
               bbox_to_anchor=(0.5, -0.10), handlelength=1.6,
               columnspacing=1.0, labelspacing=0.3)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_klevels"))


if __name__ == "__main__":
    main()
