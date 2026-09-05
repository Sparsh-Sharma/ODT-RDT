#!/usr/bin/env python3
"""The money plot for the rapid-strain follow-up: rms log-residual of the
exact-RDT-distorted vK fit (the scale-resolved ODT-vs-RDT distance) against
total strain, slow vs rapid.  Costs from the fixed-binary fit tables
(rdt_family_fit_table_S1/S20.txt); rms = sqrt(2*cost/N), N = 2*40 residuals.
"""
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
N_RES = 80.0


def read_costs(tag):
    es, costs = [], []
    for line in open(os.path.join(HERE, f"rdt_family_fit_table_{tag}.txt")):
        m = re.match(r"\s*([\d.]+)\s+[\d.]+\s+[\d.eE+-]+\s+([\d.]+)\s*$", line)
        if m:
            es.append(float(m.group(1)))
            costs.append(float(m.group(2)))
    return np.array(es), np.array(costs)


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(5.6, 4.0))
for tag, lab, col, mk in (("S1", r"$Sk/\varepsilon \approx 0.4$ (slow)", "C0", "o"),
                          ("S20", r"$Sk/\varepsilon \approx 8$ (rapid)", "C3", "s")):
    es, costs = read_costs(tag)
    rms = 100.0 * np.sqrt(2.0 * costs / N_RES)
    ax.plot(es, rms, mk + "-", color=col, ms=5, lw=1.4, label=lab)
    print(tag, [f"{r:.1f}" for r in rms])
ax.set_xlabel("total strain $e$")
ax.set_ylabel("rms log-residual of exact-RDT fit  [%]")
ax.set_title("scale-resolved ODT-vs-RDT distance\n(1023/1024 realizations, medians)")
ax.legend()
ax.set_ylim(bottom=0)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, f"fig_rdt_distance.{ext}"), dpi=200)
print("saved fig_rdt_distance.png/.pdf")
