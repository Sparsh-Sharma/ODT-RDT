#!/usr/bin/env python3
"""Manuscript figure (Part III validation subsection): four-way absolute
predictions vs measured far-field spectra, 2x2 canon panel:
(a) Paterson-Amiet NACA0012 40 m/s, (b) 90 m/s,
(c) Bampanis 2022 flat plate 32 m/s, (d) Narayanan 2015 flat plate 60 m/s.
Reuses fig_expval.run_case; writes fig_expval_paper.{pdf,png}."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..",
                                                 "closure_bound")))
from cases import CASES  # noqa: E402
from fig_expval import run_case  # noqa: E402
from figstyle_jfm import FULL, panel, plt, save  # noqa: E402

PANELS = ["PA40", "PA90", "BA22_32", "NA15_60"]
TITLES = [r"NACA0012, $U=40$ m s$^{-1}$ (Paterson--Amiet)",
          r"NACA0012, $U=90$ m s$^{-1}$ (Paterson--Amiet)",
          r"flat plate, $U=32$ m s$^{-1}$ (Bampanis $et\,al.$)",
          r"flat plate, $U=60$ m s$^{-1}$ (Narayanan $et\,al.$)"]


def main():
    fig, axs = plt.subplots(2, 2, figsize=(FULL, 4.6))
    for j, (key, ttl) in enumerate(zip(PANELS, TITLES)):
        ax = axs.flat[j]
        run_case(key, ax)
        ax.set_title(ttl, fontsize=7, pad=2.5)
        panel(ax, "abcd"[j], y=0.12)
        leg = ax.get_legend()
        if leg:
            leg.remove()
    for ax in axs[0]:
        ax.set_xlabel("")
    for ax in axs[:, 1]:
        ax.set_ylabel("")
    h, lab = axs.flat[0].get_legend_handles_labels()
    fig.legend(h, lab, ncol=5, loc="lower center", fontsize=6.5,
               bbox_to_anchor=(0.5, -0.06), columnspacing=1.0)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_expval_paper"))


if __name__ == "__main__":
    main()
