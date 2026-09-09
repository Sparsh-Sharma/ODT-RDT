#!/usr/bin/env python3
"""Manuscript figure for sec:envelope (main_v3.tex): rms log-residual of
the exact-RDT-distorted vK fit against total strain, slow vs rapid.
Costs from rdt_family_fit_table_S1/S20.txt (mirrored from LEN_Extension
post/acoustics/results_test3, fixed-binary ensembles 2026-09-05);
rms = sqrt(2*cost/N), N = 2*40 residuals.  Output at 0.62\\textwidth."""
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from figstyle_jfm import FULL, plt, save  # noqa: E402

N_RES = 80.0


def read_costs(tag):
    es, costs = [], []
    for line in open(os.path.join(HERE, f"rdt_family_fit_table_{tag}.txt")):
        m = re.match(r"\s*([\d.]+)\s+[\d.]+\s+[\d.eE+-]+\s+([\d.]+)\s*$",
                     line)
        if m:
            es.append(float(m.group(1)))
            costs.append(float(m.group(2)))
    return np.array(es), np.array(costs)


def main():
    fig, ax = plt.subplots(figsize=(0.62 * FULL, 2.3))
    # monochrome: slow -> solid + filled circle, rapid -> dashed + filled
    # square (rapidity by line style + marker).
    for tag, lab, ls, mk in (
            ("S1", r"$Sk_t/\varepsilon \approx 0.4$ (slow)", "-", "o"),
            ("S20", r"$Sk_t/\varepsilon \approx 8$ (rapid)",
             (0, (5, 2)), "s")):
        es, costs = read_costs(tag)
        rms = 100.0 * np.sqrt(2.0 * costs / N_RES)
        ax.plot(es, rms, color="k", ls=ls, lw=1.1, marker=mk, ms=3.4,
                mfc="k", mew=0.6, label=lab)
        print(tag, [f"{r:.1f}" for r in rms])
    ax.set_xlabel("total strain $e$")
    ax.set_ylabel("rms log-residual of exact-RDT fit [%]")
    ax.legend(loc="lower right")
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_rdt_distance"))


if __name__ == "__main__":
    main()
