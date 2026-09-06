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
from figstyle_jfm import FULL, plt, save, strip  # noqa: E402

NB = os.path.join(HERE, "newbands")
CASES = [
    ("standard ODT", "bands_homogeneousStrain2.npz", "k", "o"),
    ("1 level (thirds)", "bands_homogeneousStrain2K1.npz", "#2a78d6", "s"),
    ("2 levels", "bands_homogeneousStrain2K2.npz", "#1baf7a", "^"),
    ("3 levels", "bands_homogeneousStrain2K3.npz", "#eb6834", "d"),
    ("depth tied to eddy size", "bands_homogeneousStrain2KS.npz",
     "#8e44ad", "v"),
    ("2 levels, iterated 3x", "bands_homogeneousStrain2K2I3.npz",
     "#b8860b", "P"),
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
    for lab, fn, col, mk in CASES:
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
            ax.errorbar(es, m, yerr=[m - lo, hi - m], fmt=mk + "-",
                        color=col, capsize=1.8, capthick=0.6,
                        elinewidth=0.6, ms=3, lw=0.9, label=lab)
            ax.set_ylabel(ylab)
    axs[2].axhline(1.0, color="k", lw=0.6, ls=":")
    for j, ax in enumerate(axs):
        strip(ax)
        if j >= 2:
            ax.set_xlabel("total strain $e$")
        ax.text(0.03, 0.97, "(" + "abcd"[j] + ")",
                transform=ax.transAxes, va="top")
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, -0.10), handlelength=1.6,
               columnspacing=1.2, labelspacing=0.3)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_klevels"))


if __name__ == "__main__":
    main()
