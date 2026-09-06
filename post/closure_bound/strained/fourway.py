#!/usr/bin/env python3
"""Four-way comparison (main_v3.tex sec:limits capstone): the
strain-induced upwash anisotropy per wavenumber band, Delta b_22(k2),
for exact linear RDT, DNS, the strain-coupled model (ISO baseline), and
the strain-coupled model with the adaptive-depth sub-scale kernels on
(the scale-local-relaxation variant) — at Sk/eps ~ 0.8 (S=2 ensembles)
and Sk/eps ~ 16 (S=40), all at e=1, on the observables of threeway.py.

ODT sides: fresh same-binary ensembles fw_{S2,S40}{I,K} (caro
2026-09-08, 1024 rlz each, gateA_S1 precursor protocol).  The fresh ISO
baseline is cross-checked against the archived s2_iso ensemble.
Usage: python fourway.py <dir with fw_*_ensemble.npz>
Writes fourway.npz + fig_fourway.{pdf,png} (canon).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, os.pardir))
sys.path.insert(0, os.path.join(HERE, os.pardir, "odt_alloc"))
os.environ.setdefault("ALLOC_SLOW", "S2")
import threeway as T  # noqa: E402  (D, EDGES, XC, observables, dns_like)
from alloc_tests import centroid, line_spectra  # noqa: E402
from figstyle_jfm import COL, FULL, GREEN, panel, plt, save  # noqa: E402


def odt_ens(npz_path, di, f):
    d = np.load(npz_path)
    k0, (q11, _), (q22, _), (q33, _) = line_spectra(
        d["lines"][0], d["Ldump"][0], all_components=True)
    k1, (p11, _), (p22, _), (p33, _) = line_spectra(
        d["lines"][di], d["Ldump"][di], all_components=True)
    kref = centroid(k0, q11 + q22 + q33)
    return T.observables(k0, [q11, q22, q33], k1, [p11, p22, p33], f,
                         12.0 * kref)


def main():
    fwdir = sys.argv[1] if len(sys.argv) > 1 else HERE
    f = np.exp(0.5)          # e = 1
    res = {}
    for rat, tag in (("0.8", "S2"), ("16", "S40")):
        res[(rat, "RDT")] = T.dns_like(f"chk_r{rat}_s*_e0_rdt.npz",
                                       f"chk_r{rat}_s*_e1_rdt.npz", f)
        res[(rat, "DNS")] = T.dns_like(f"chk_r{rat}_s*_e0.npz",
                                       f"chk_r{rat}_s*_e1.npz", f)
        for suff, lab in (("I", "ODT"), ("K", "ODTK")):
            p = os.path.join(fwdir, f"fw_{tag.lower()}{suff.lower()}"
                                    "_ensemble.npz")
            res[(rat, lab)] = odt_ens(p, 4, f)
    # cross-check: fresh ISO baseline vs archived s2_iso ensemble
    arch = os.path.join(HERE, os.pardir, "odt_alloc",
                        "s2_iso_ensemble.npz")
    if os.path.exists(arch):
        old = odt_ens(arch, 4, f)
        dd = np.nanmax(np.abs(old["db22"] - res[("0.8", "ODT")]["db22"]))
        print(f"ISO baseline fresh-vs-archived max |d db22| = {dd:.4f}")

    np.savez(os.path.join(HERE, "fourway.npz"), edges=T.EDGES,
             **{f"r{rat}_{lab}_{q}": res[(rat, lab)][q]
                for (rat, lab) in res for q in ("db22", "db11", "db33")})

    CC = {"RDT": COL["rdt"], "DNS": COL["dns"], "ODT": COL["odt"],
          "ODTK": GREEN}
    LBL = {"RDT": "exact linear RDT", "DNS": r"DNS $128^3$",
           "ODT": "strain-coupled ODT",
           "ODTK": "ODT $+$ scale-local relaxation"}
    MK = {"RDT": "D", "DNS": "o", "ODT": "s", "ODTK": "^"}
    fig, axs = plt.subplots(1, 2, figsize=(FULL, 2.3), sharey=True)
    for ax, rat, ttl in ((axs[0], "0.8", r"$Sk_t/\varepsilon=0.8$"),
                         (axs[1], "16", r"$Sk_t/\varepsilon=16$")):
        for lab in ("RDT", "DNS", "ODT", "ODTK"):
            ax.plot(T.XC, res[(rat, lab)]["db22"], "-", color=CC[lab],
                    marker=MK[lab], ms=2.8,
                    label=LBL[lab] if ax is axs[0] else None)
        ax.axhline(0, color="0.6", lw=0.6, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel(r"$\kappa_2(e)/\kappa_c(0)$")
        ax.text(0.5, 0.96, ttl, transform=ax.transAxes, ha="center",
                va="top", fontsize=8)
    axs[0].set_ylabel(r"$\Delta b_{22}(\kappa_2)$")
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, ncol=4, loc="lower center",
               bbox_to_anchor=(0.5, -0.14), columnspacing=1.0)
    for j, ax in enumerate(axs):
        panel(ax, "ab"[j], y=0.14)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_fourway"))
    for rat in ("0.8", "16"):
        print(f"-- Sk/eps={rat}: db22 bands")
        for lab in ("RDT", "DNS", "ODT", "ODTK"):
            v = res[(rat, lab)]["db22"]
            print(f"  {lab:5s} " + " ".join(f"{x:+.3f}" for x in v))


if __name__ == "__main__":
    main()
