#!/usr/bin/env python3
"""anisoRejectFac sweep summary (baseline + fac 0.9/0.8/0.7/0.5, 1024 rlz each,
paired seeds).  Robust stats: median over realizations, bootstrap CI of median.
Writes facSweep_table.txt and fig_facSweep.{png,pdf} next to this script.
"""
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = [("baseline", "bands_homogeneousStrain2.npz", None, "k"),
         ("fac 0.9", "bands_homogeneousStrain2A.npz", 45, "C0"),
         ("fac 0.8", "bands_homogeneousStrain2A80.npz", 66, "C2"),
         ("fac 0.7", "bands_homogeneousStrain2A70.npz", 83, "C1"),
         ("fac 0.5", "bands_homogeneousStrain2A50.npz", 98, "C3"),
         ("S 0.9", "bands_homogeneousStrain2AS90.npz", 42, "C4"),   # scale-conditioned, l*=0.05
         ("S 0.5", "bands_homogeneousStrain2AS50.npz", 97, "C5")]   # rej% of gated (sub-l*) candidates
NBOOT = 4000


def med_ci(x, rng):
    m = np.median(x)
    idx = rng.integers(0, x.size, (NBOOT, x.size))
    mb = np.median(x[idx], axis=1)
    return m, np.percentile(mb, 2.5), np.percentile(mb, 97.5)


def main():
    rng = np.random.default_rng(0)
    data = {}
    for lab, fn, rej, col in CASES:
        d = np.load(os.path.join(HERE, fn))
        data[lab] = {"Alo": d["E2lo"] / d["Eplo"], "Ahi": d["E2hi"] / d["Ephi"],
                     "T": (d["E2hi"] / d["Ephi"]) / (d["E2lo"] / d["Eplo"]),
                     "u2": d["u2frac"], "strains": d["strains"], "rej": rej,
                     "col": col}

    strains = data["baseline"]["strains"]
    lines = ["anisoRejectFac sweep, 1024 paired-seed realizations per case",
             "median over realizations [95% bootstrap CI of median]", ""]
    for j, e in enumerate(strains):
        if e == 0:
            continue
        lines.append(f"--- e = {e:.1f} ---")
        lines.append(f"{'case':>9} {'rej%':>5} {'A_low':>19} {'A_high':>19}"
                     f" {'T=A_hi/A_lo':>19} {'u2/2kt':>7}")
        for lab, *_ in CASES:
            d = data[lab]
            alo = med_ci(d["Alo"][:, j], rng)
            ahi = med_ci(d["Ahi"][:, j], rng)
            t = med_ci(d["T"][:, j], rng)
            rej = f"{d['rej']:d}" if d["rej"] else "0"
            lines.append(
                f"{lab:>9} {rej:>5} "
                f"{alo[0]:5.2f} [{alo[1]:4.2f},{alo[2]:4.2f}] "
                f"{ahi[0]:5.2f} [{ahi[1]:4.2f},{ahi[2]:4.2f}] "
                f"{t[0]:5.2f} [{t[1]:4.2f},{t[2]:4.2f}] "
                f"{np.median(d['u2'][:, j]):7.3f}")
        lines.append("")
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(HERE, "facSweep_table.txt"), "w").write(txt + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    es = strains
    for lab, *_ in CASES:
        d = data[lab]
        for ax, key in zip(axes, ("Alo", "Ahi", "T")):
            m = np.array([med_ci(d[key][:, j], rng) for j in range(len(es))])
            ax.errorbar(es, m[:, 0], yerr=[m[:, 0] - m[:, 1], m[:, 2] - m[:, 0]],
                        fmt="o-", color=d["col"], ms=4, capsize=2, lw=1.2,
                        label=lab + (f" ({d['rej']}% rej)" if d["rej"] else ""))
    for ax, ttl in zip(axes, (r"$A_{\rm low}$ (k2 30-100)",
                              r"$A_{\rm high}$ (k2 300-800)",
                              r"transmission $A_{\rm high}/A_{\rm low}$")):
        ax.axhline(1.0, color="k", lw=0.8, ls=":")
        ax.set_xlabel("total strain e")
        ax.set_title(ttl)
    axes[0].set_ylabel("median over realizations, 95% CI")
    axes[0].legend(fontsize=8)
    fig.suptitle("Option A (anisotropy-gated eddy acceptance): threshold sweep, 1024 rlz/case")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"fig_facSweep.{ext}"), dpi=200)
    print("figure saved: fig_facSweep.png/.pdf")


if __name__ == "__main__":
    main()
