#!/usr/bin/env python3
"""Five-way extension of the four-way comparison: add the size-capped
relaxation-clock variant (fw_S1R/S2R/S40R, rate 300 matched to the
strained-phase eddy rate, Lmax 0.05, depth 6) to the Delta b_22(k2)
axes at all three rapidities.  Extends fourway.npz with r*_ODTR rows
and renders fig_fourway5.{pdf,png}.
Usage: python fourway_rcs.py <dir with fw_s*r_ensemble.npz>
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, os.pardir))
sys.path.insert(0, os.path.join(HERE, os.pardir, "odt_alloc"))
os.environ.setdefault("ALLOC_SLOW", "S2")
import threeway as T  # noqa: E402
from fourway import odt_ens  # noqa: E402
from figstyle_jfm import COL, FULL, GREEN, RED, panel, plt, save  # noqa: E402


def main():
    fwdir = sys.argv[1]
    f = np.exp(0.5)
    d = dict(np.load(os.path.join(HERE, "fourway.npz")))
    for rat, tag in (("0.4", "s1"), ("0.8", "s2"), ("16", "s40")):
        o = odt_ens(os.path.join(fwdir, f"fw_{tag}r_ensemble.npz"), 4, f)
        for q in ("db22", "db11", "db33"):
            d[f"r{rat}_ODTR_{q}"] = o[q]
    np.savez(os.path.join(HERE, "fourway.npz"), **d)

    CC = {"RDT": COL["rdt"], "DNS": COL["dns"], "ODT": COL["odt"],
          "ODTK": GREEN, "ODTR": RED}
    LBL = {"RDT": "exact linear RDT", "DNS": r"DNS $128^3$",
           "ODT": "strain-coupled ODT",
           "ODTK": "ODT $+$ scale-local relaxation",
           "ODTR": "ODT $+$ capped relaxation clock"}
    MK = {"RDT": "D", "DNS": "o", "ODT": "s", "ODTK": "^", "ODTR": "v"}
    fig, axs = plt.subplots(1, 3, figsize=(FULL, 2.1), sharey=True)
    for ax, rat, ttl in ((axs[0], "0.4", r"$Sk_t/\varepsilon=0.4$"),
                         (axs[1], "0.8", r"$Sk_t/\varepsilon=0.8$"),
                         (axs[2], "16", r"$Sk_t/\varepsilon=16$")):
        for lab in ("RDT", "DNS", "ODT", "ODTK", "ODTR"):
            key = f"r{rat}_{lab}_db22"
            if key not in d:
                if lab == "RDT":
                    key = "r0.8_RDT_db22"
                else:
                    continue
            ax.plot(T.XC, d[key], "-", color=CC[lab], marker=MK[lab],
                    ms=2.4, label=LBL[lab] if ax is axs[1] else None)
        ax.axhline(0, color="0.6", lw=0.6, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel(r"$\kappa_2(e)/\kappa_c(0)$")
        ax.text(0.5, 0.96, ttl, transform=ax.transAxes, ha="center",
                va="top", fontsize=7.5)
    axs[0].set_ylabel(r"$\Delta b_{22}(\kappa_2)$")
    axs[0].text(0.5, 0.85, "(no DNS at 0.4)", transform=axs[0].transAxes,
                ha="center", fontsize=6.5, color="0.4")
    h, lab = axs[1].get_legend_handles_labels()
    fig.legend(h, lab, ncol=5, loc="lower center",
               bbox_to_anchor=(0.5, -0.16), columnspacing=0.8,
               fontsize=6.3)
    for j, ax in enumerate(axs):
        panel(ax, "abc"[j], y=0.14)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_fourway5"))
    for rat in ("0.4", "0.8", "16"):
        print(f"-- Sk/eps={rat}")
        for lab in ("ODT", "ODTK", "ODTR"):
            k = f"r{rat}_{lab}_db22"
            if k in d:
                print(f"  {lab:5s} "
                      + " ".join(f"{x:+.3f}" for x in d[k]))


if __name__ == "__main__":
    main()
