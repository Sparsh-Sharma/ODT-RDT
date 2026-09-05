#!/usr/bin/env python3
"""Gate A refit on the 1024-realization CARO baseline ensemble (median
spectra), then Gate B Delta-SPL.

Targets: per-k2 MEDIAN ensemble spectra from dump_spectra.py
(spectra_homogeneousStrain2.npz).  For each strain e the axisymmetric
vK-stretch family (A0, c0, L2, Lperp) is fitted to (Eperp, E2); Gate B then
compares the strained fits against
  (a) the fitted e=0 spectrum      -> total Delta-SPL (incl. amplification),
  (b) matched_isotropic(mode=...)  -> shape-only Delta-SPL.

Decisions taken here (were DECISION NEEDED):
  - h(mu): P2-Legendre stub kept (paper-1 kernel source not in repo).
  - C-tensor: rank-1 projected-axis dyad kept.
  - ell_y: Amiet reciprocal correlation length (gateB CALL 1) — unchanged.
  - Baseline: BOTH the measured e=0 fit (total) and matched_isotropic
    (shape-only) are reported; the e=0 fit is the physical reference.

Writes gateA_fit_table.txt and fig_gateA_1024.{png,pdf} next to this script.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import odt_io                                          # noqa: E402
import fit_family as ff                                # noqa: E402
import gateB_delta_spl as gb                           # noqa: E402
from axisym_family import ke_from_length               # noqa: E402

NPZ = os.path.join(HERE, "spectra_homogeneousStrain2.npz")
NBIN = 40
KMAX = 1200.0


def fit_one(k2, Eperp, E2, kmin):
    sel = (k2 >= kmin) & (k2 <= KMAX)
    k2b, (Epb, E2b) = odt_io.log_bin(k2[sel], [Eperp[sel], E2[sel]], NBIN)
    res = ff.fit_family(k2b, Epb, E2b)
    return res, k2b, Epb, E2b


def main():
    d = np.load(NPZ)
    strains = d["strains"]
    lines = [f"Gate A refit, median ensemble spectra, {int(d['nrlz'])} rlz "
             f"(homogeneousStrain2)",
             f"band: k2 in [3*dk, {KMAX:.0f}], log-binned to {NBIN} pts; "
             f"h(mu)=P2 stub, rank-1 C-tensor",
             "",
             f"{'e':>5} {'A0':>11} {'c0':>8} {'L2':>9} {'Lperp':>9} "
             f"{'Lperp/L2':>9} {'cost':>9} {'E1/E3':>7}"]
    fits, targets = [], []
    for j, e in enumerate(strains):
        k2 = d[f"k2_e{j}"]
        E1, E2, E3 = d[f"med_E1_e{j}"], d[f"med_E2_e{j}"], d[f"med_E3_e{j}"]
        Eperp = 0.5 * (E1 + E3)
        kmin = 3.0 * k2[0]
        res, k2b, Epb, E2b = fit_one(k2, Eperp, E2, kmin)
        p = res.params
        m = (k2 >= kmin) & (k2 <= KMAX)
        r13 = np.median(E1[m] / np.clip(E3[m], 1e-300, None))
        fits.append(p)
        targets.append((k2b, Epb, E2b))
        lines.append(f"{e:5.1f} {p.A0:11.4e} {p.c0:8.3f} {p.L2:9.3f} "
                     f"{p.Lperp:9.3f} {p.Lperp/p.L2:9.3f} {res.cost:9.3f} "
                     f"{r13:7.3f}")

    # ---- Gate B ----
    p_iso = fits[0]
    ke0 = max(p_iso.L2, p_iso.Lperp)
    kx = np.geomspace(0.2, 40.0, 13) * ke0
    lines += ["", "Gate B Delta-SPL [dB] vs K_x/k_e0 "
              "(total = vs fitted e=0; shape = vs matched_isotropic 'perp'):"]
    hdr = f"{'K_x/k_e0':>9}"
    for e in strains[1:]:
        hdr += f"  tot e={e:<4.1f} shp e={e:<4.1f}"
    lines.append(hdr)
    curves = {}
    for j, e in enumerate(strains[1:], start=1):
        p = fits[j]
        tot = np.array([10*np.log10(gb.noise_kernel(k, p) /
                                    gb.noise_kernel(k, p_iso)) for k in kx])
        p_m = gb.matched_isotropic(p, mode="perp")
        shp = np.array([10*np.log10(gb.noise_kernel(k, p) /
                                    gb.noise_kernel(k, p_m)) for k in kx])
        curves[e] = (tot, shp)
    for i, k in enumerate(kx):
        row = f"{k/ke0:9.2f}"
        for e in strains[1:]:
            row += f"  {curves[e][0][i]:9.2f} {curves[e][1][i]:9.2f}"
        lines.append(row)

    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(HERE, "gateA_fit_table.txt"), "w").write(txt + "\n")

    # ---- figure ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.0))
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(strains)))
    for j, (e, p, (k2b, Epb, E2b)) in enumerate(zip(strains, fits, targets)):
        E1m, E2m = ff.line_spectra(p, k2b)
        axs[0].loglog(k2b, E2b, "o", ms=2.5, color=cols[j])
        axs[0].loglog(k2b, E2m, "-", lw=1.1, color=cols[j],
                      label=f"e={e:.1f}")
        axs[1].loglog(k2b, Epb, "o", ms=2.5, color=cols[j])
        axs[1].loglog(k2b, E1m, "-", lw=1.1, color=cols[j])
    axs[0].set_title("$E_2(k_2)$ (upwash), median + fit")
    axs[1].set_title(r"$E_\perp(k_2)$, median + fit")
    for ax in axs[:2]:
        ax.set_xlabel("$k_2$")
    axs[0].legend(fontsize=7)
    for j, e in enumerate(strains[1:], start=1):
        axs[2].semilogx(kx/ke0, curves[e][0], "-", color=cols[j],
                        label=f"total, e={e:.1f}")
        axs[2].semilogx(kx/ke0, curves[e][1], "--", color=cols[j], lw=1.0)
    axs[2].axhline(0, color="k", lw=0.7, ls=":")
    axs[2].set_xlabel("$K_x/k_{e0}$")
    axs[2].set_ylabel(r"$\Delta$SPL [dB]")
    axs[2].set_title("Gate B (solid: vs e=0 fit; dashed: shape-only)")
    axs[2].legend(fontsize=7)
    fig.suptitle("Gate A/B on 1024-rlz median spectra (homogeneousStrain2)",
                 fontsize=11)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"fig_gateA_1024.{ext}"), dpi=200)
    print("saved fig_gateA_1024.png/.pdf")


if __name__ == "__main__":
    main()
