#!/usr/bin/env python3
"""Alan's last-eddy-size clock variant (2026-09-11) against the existing
clock family on the transmission diagnostic.  (a) fine-scale-band imbalance
A_hi(e); (b) upwash energy fraction u2^2/2kt(e).  Cases: baseline, uncapped
clock (RC1, distribution sizes), size-capped clock (RCS1), last-eddy-size
clock uncapped (RCE1) and capped (RCES1).  All 1024 paired-seed realisations,
one binary per campaign pair.  Email figure (colour)."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NB = os.path.join(HERE, "newbands")

CASES = [
    ("SC-ODT (baseline)", "bands_homogeneousStrain2.npz", "o", "-", "k"),
    ("clock, uncapped (RC1)", "bands_homogeneousStrain2RC1.npz", "s", "--",
     "tab:blue"),
    ("clock, size-capped (RCS1)", "bands_homogeneousStrain2RCS1.npz", "^",
     "-.", "tab:green"),
    ("last-eddy clock, uncapped (RCE1)", "bands_homogeneousStrain2RCE1.npz",
     "D", (0, (1, 1)), "tab:red"),
    ("last-eddy clock, capped (RCES1)", "bands_homogeneousStrain2RCES1.npz",
     "v", (0, (4, 1, 1, 1)), "tab:orange"),
]
NBOOT = 2000


def med_ci(x, rng):
    m = np.median(x)
    idx = rng.integers(0, x.size, (NBOOT, x.size))
    mb = np.median(x[idx], axis=1)
    return m, np.percentile(mb, 2.5), np.percentile(mb, 97.5)


def main():
    rng = np.random.default_rng(0)
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 3.2), sharex=True)
    for lab, fn, mk, ls, col in CASES:
        d = np.load(os.path.join(NB, fn))
        es = d["strains"]
        Ahi = d["E2hi"] / d["Ephi"]
        m, lo, hi = np.array([med_ci(Ahi[:, j], rng)
                              for j in range(len(es))]).T
        a.errorbar(es, m, yerr=[m - lo, hi - m], fmt=mk, ls=ls, color=col,
                   capsize=2, ms=4, lw=1.2, label=lab)
        u2 = np.median(d["u2frac"], axis=0)
        b.plot(es, u2, marker=mk, ls=ls, color=col, ms=4, lw=1.2)
    a.axhline(1.0, color="0.7", lw=0.6, zorder=0)
    a.set_ylabel(r"$\mathcal{A}_{\rm hi}$ (fine-scale band imbalance)")
    a.set_xlabel(r"total strain $e$")
    b.set_ylabel(r"$u_2^2/2k_t$")
    b.set_xlabel(r"total strain $e$")
    h, lab = a.get_legend_handles_labels()
    fig.legend(h, lab, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, -0.12), fontsize=8)
    a.set_title("(a) fine scales: hold vs fade", fontsize=9)
    b.set_title("(b) one-point price", fontsize=9)
    fig.tight_layout(pad=0.5)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"fig_clock_lastEddy.{ext}"),
                    dpi=180, bbox_inches="tight")
    print("wrote fig_clock_lastEddy.png/pdf")


if __name__ == "__main__":
    main()
