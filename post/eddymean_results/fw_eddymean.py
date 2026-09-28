#!/usr/bin/env python3
"""Eddy-on-mean re-validation (T2, decisive spectral test): Delta b_22(k2) at
e=1 for exact RDT, 128^3 DNS, baseline SC-ODT, gated eddy-on-mean (EM), and
ungated eddy-on-mean (EMU, mechanism ceiling), at Sk/eps~0.8 (S2) and ~16
(S40). Archived adaptive-depth (ODTK) overlaid from fourway.npz.
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
from alloc_tests import centroid, line_spectra  # noqa: E402
from figstyle_jfm import FULL, panel, plt, save  # noqa: E402


def odt_ens(npz_path, di, f):
    d = np.load(npz_path)
    k0, (q11, _), (q22, _), (q33, _) = line_spectra(
        d["lines"][0], d["Ldump"][0], all_components=True)
    k1, (p11, _), (p22, _), (p33, _) = line_spectra(
        d["lines"][di], d["Ldump"][di], all_components=True)
    return T.observables(k0, [q11, q22, q33], k1, [p11, p22, p33], f,
                         12.0 * centroid(k0, q11 + q22 + q33))


VARIANTS = [("I", "baseline SC-ODT", "s", "#2a78d6"),
            ("EM", "eddy-on-mean, gated", "^", "#1baf7a"),
            ("EMU", "eddy-on-mean, ungated", "v", "#c22f2f")]


def main():
    f = np.exp(0.5)          # e = 1
    res = {}
    for rat, tag in (("0.8", "s2"), ("16", "s40")):
        res[(rat, "RDT")] = T.dns_like(f"chk_r{rat}_s*_e0_rdt.npz",
                                       f"chk_r{rat}_s*_e1_rdt.npz", f)
        res[(rat, "DNS")] = T.dns_like(f"chk_r{rat}_s*_e0.npz",
                                       f"chk_r{rat}_s*_e1.npz", f)
        for suff, _, _, _ in VARIANTS:
            p = os.path.join(HERE, f"fw_{tag}{suff.lower()}_ensemble.npz")
            res[(rat, suff)] = odt_ens(p, 4, f)

    arch = np.load(os.path.join(HERE, "fourway.npz"))
    for rat in ("0.8", "16"):
        k = f"r{rat}_ODTK_db22"
        if k in arch.files:
            res[(rat, "ODTK")] = {"db22": arch[k]}

    print("bands k2(e)/kc(0): " +
          " ".join(f"[{a:.2f},{b:.2f})" for a, b in zip(T.EDGES[:-1],
                                                        T.EDGES[1:])))
    order = ["RDT", "DNS", "I", "EM", "EMU", "ODTK"]
    for rat in ("0.8", "16"):
        print(f"\n===== Sk/eps ~ {rat} =====  db22(k2) per band, e=1")
        for s in order:
            if (rat, s) in res:
                v = res[(rat, s)]["db22"]
                print(f"  {s:5s} " + " ".join(f"{x:+6.3f}" for x in v))

    fig, axs = plt.subplots(1, 2, figsize=(FULL, 2.4), sharey=True)
    for ax, rat, ttl in ((axs[0], "0.8", r"$Sk_t/\varepsilon=0.8$"),
                         (axs[1], "16", r"$Sk_t/\varepsilon=16$")):
        ax.plot(T.XC, res[(rat, "RDT")]["db22"], "-", color="#8e44ad",
                marker="D", ms=2.6, label="exact RDT")
        ax.plot(T.XC, res[(rat, "DNS")]["db22"], "-", color="k",
                marker="o", ms=2.8, label=r"DNS $128^3$")
        for suff, lbl, mk, cc in VARIANTS:
            ax.plot(T.XC, res[(rat, suff)]["db22"], "-", color=cc,
                    marker=mk, ms=2.8, label=lbl)
        if (rat, "ODTK") in res:
            ax.plot(T.XC, res[(rat, "ODTK")]["db22"], ":", color="0.5",
                    lw=0.9, label="adaptive-depth (arch.)")
        ax.axhline(0, color="0.6", lw=0.6, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel(r"$\kappa_2(e)/\kappa_c(0)$")
        ax.text(0.5, 0.05, ttl, transform=ax.transAxes, ha="center",
                fontsize=8)
    axs[0].set_ylabel(r"$\Delta b_{22}(\kappa_2)$")
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, -0.20), columnspacing=1.0, fontsize=7)
    for j, ax in enumerate(axs):
        panel(ax, "ab"[j], y=0.96)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_fw_eddymean"))


if __name__ == "__main__":
    main()
