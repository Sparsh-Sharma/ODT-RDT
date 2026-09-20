#!/usr/bin/env python3
"""Summary ("hero") figure for the SC-ODT programme.

One left-to-right arc, all from real ensembles already in the repo:
  (a) distortion  -- measured upwash line spectra at e = 0, 1, 2
                     (spectra_gateA_S1.npz, ~1024-realisation medians);
  (b) validation  -- strain-induced upwash anisotropy per wavenumber,
                     SC-ODT against the 128^3 strained-box DNS and exact
                     projected RDT at e = 1 (strained/threeway.npz);
  (c) consequence -- signed leading-edge-noise prediction error against
                     total distortion strain (lenoise/fig_frozen_cost).

Two renders from the same data:
  python fig_hero.py            -> monochrome, JFM paper canon
  python fig_hero.py --color    -> colour, for proposals

No hype; panel titles state what is shown, axes state the quantity.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "lenoise"))
sys.path.insert(0, os.path.join(HERE, "strained"))

import figstyle_jfm as S            # noqa: E402
from figstyle_jfm import FULL, panel, plt, save, COL  # noqa: E402

COLOR = "--color" in sys.argv

# ---------------------------------------------------------------- palettes
if COLOR:
    RAMP = {0: "0.35", 1: S.ORANGE, 2: S.RED}     # e=0 grey (reference)
    SYS = {"DNS": COL["dns"], "RDT": COL["rdt"], "ISO": COL["odt"]}
    CSER = {"vk": "k", "std": S.BLUE, "fix": S.GREEN}
    FILL0 = S.BLUE
else:
    RAMP = {0: "0.62", 1: "0.35", 2: "k"}
    SYS = {"DNS": "k", "RDT": "k", "ISO": "k"}
    CSER = {"vk": "k", "std": "0.45", "fix": "0.45"}
    FILL0 = "0.85"

MK_SYS = {"RDT": "D", "DNS": "o", "ISO": "s"}
LBL_SYS = {"RDT": "exact RDT", "DNS": r"DNS $128^3$", "ISO": "SC-ODT"}


# --------------------------------------------------------------- panel (a)
def panel_distortion(ax):
    d = np.load(os.path.join(HERE, "lenoise", "spectra_gateA_S1.npz"),
                allow_pickle=True)
    idx = {0: 0, 1: 2, 2: 4}          # dump index = 2 * e  (e = 0, 1, 2)
    lab = {0: "isotropic (frozen inflow)", 1: "$e=1$", 2: "strained, $e=2$"}
    curves = {}
    for e, i in idx.items():
        k = d[f"k2_e{i}"]
        E = d[f"med_E2_e{i}"]
        m = (k > 0) & (E > 0)
        curves[e] = (k[m], E[m])
    # shade what strain does: between the isotropic and the e=2 spectra
    k0, E0 = curves[0]
    k2, E2 = curves[2]
    kk = np.geomspace(max(k0.min(), k2.min()), min(k0.max(), k2.max()), 400)
    f0 = np.interp(kk, k0, E0)
    f2 = np.interp(kk, k2, E2)
    ax.fill_between(kk, f0, f2, color=FILL0, alpha=0.16, lw=0, zorder=0)
    for e in (0, 1, 2):
        k, E = curves[e]
        lw = 1.6 if e in (0, 2) else 1.1
        ls = "-" if e in (0, 2) else (0, (5, 1.8))
        ax.plot(k, E, color=RAMP[e], lw=lw, ls=ls, label=lab[e],
                zorder=3 if e else 4)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"wavenumber $\kappa_2$")
    ax.set_ylabel(r"upwash spectrum $E_{2}(\kappa_2)$")
    ax.set_title("incoming turbulence, distorted by strain", pad=4)
    ax.legend(loc="lower left", fontsize=6.6)
    # arrow: energy-containing scales shift to higher wavenumber
    ax.annotate("", xy=(k2[np.argmax(E2)] * 1.0, E2.max()),
                xytext=(k0[np.argmax(E0)] * 1.0, E0.max()),
                arrowprops=dict(arrowstyle="-|>", color="0.4", lw=0.8,
                                shrinkA=2, shrinkB=2), zorder=5)


# --------------------------------------------------------------- panel (b)
def panel_validation(ax):
    D = np.load(os.path.join(HERE, "strained", "threeway.npz"),
                allow_pickle=True)
    E = D["edges"]
    xc = np.sqrt(E[:-1] * E[1:])
    ax.axhline(0, color="0.75", lw=0.6, ls=":")
    for s in ("DNS", "RDT", "ISO"):
        y = D[f"{s}_e1.0_db22"]
        ax.plot(xc, y, color=SYS[s], ls="-", lw=1.2 if s == "DNS" else 1.0,
                marker=MK_SYS[s], ms=3.4, mfc=SYS[s], mew=0.6,
                label=LBL_SYS[s], zorder=4 if s == "ISO" else 3)
    ax.set_xscale("log")
    ax.set_xlabel(r"wavenumber $\kappa_2(e)/\kappa_c(0)$")
    ax.set_ylabel(r"upwash anisotropy $\Delta b_{22}$")
    ax.set_title(r"validated against DNS, $e=1$", pad=4)
    ax.legend(loc="upper right", fontsize=6.6)


# --------------------------------------------------------------- panel (c)
def panel_consequence(ax):
    import fig_frozen_cost as ffc     # noqa: E402
    import fig_amiet_rapid as far     # noqa: E402
    app = ffc.app_errors()
    ax.axvspan(0.020, 0.29, color="0.94", zorder=0)
    ax.axhline(0.0, color="0.75", lw=0.6, zorder=1)
    order = ("vk", "std", "fix")
    mk = {"vk": "s", "std": "o", "fix": "^"}
    ls = {"vk": (0, (5, 2)), "std": "-", "fix": (0, (1, 1.4))}
    lab = {"vk": "frozen von Kármán", "std": "computed (SC-ODT)",
           "fix": "computed (clock-relaxed)"}
    for key in order:
        col = CSER[key]
        for rows, filled in ((ffc.PLATES, True), (ffc.NACA, False)):
            ci = {"vk": 4, "std": 5, "fix": 6}[key]
            j = ffc.JIT[key]
            x = np.array([r[1] for r in rows]) * j
            xe = np.array([[r[1] - r[2] for r in rows],
                           [r[3] - r[1] for r in rows]]) * j
            yv = [r[ci] for r in rows]
            ax.errorbar(x, yv, xerr=xe, fmt=mk[key], color=col, ms=3.0,
                        mfc=(col if filled else "white"), mec=col,
                        ecolor="0.82", elinewidth=0.5, capsize=0,
                        ls="none", zorder=4)
        ax.plot(ffc.ES_APP, app[key], color=col, ls=ls[key], marker=mk[key],
                ms=3.0, mfc=col, lw=1.2, zorder=5, label=lab[key])
        if key in far.RMS:
            eb = [e for e in ffc.ES_APP if e in far.RMS[key]]
            db = np.array([10 * np.log10(1 + far.RMS[key][e]) for e in eb])
            yb = np.array([app[key][ffc.ES_APP.index(e)] for e in eb])
            ax.fill_between(eb, yb - db, yb + db,
                            color=(col if COLOR else "0.85"),
                            alpha=0.16 if COLOR else 0.55, lw=0, zorder=2)
    ax.annotate("cost of freezing the inflow",
                xy=(1.9, app["vk"][-1]), xytext=(0.33, 5.1), fontsize=6.6,
                ha="left", arrowprops=dict(arrowstyle="-|>", color="k",
                                           lw=0.7, shrinkB=3))
    ax.set_xscale("log")
    ax.set_xlim(0.020, 3.2)
    ax.set_ylim(-4.4, 6.9)
    ax.set_xticks([0.03, 0.1, 0.3, 1, 2])
    ax.set_xticklabels(["0.03", "0.1", "0.3", "1", "2"])
    ax.set_xlabel(r"total distortion strain $e$")
    ax.set_ylabel(r"noise prediction error [dB]")
    ax.set_title("effect on predicted leading-edge noise", pad=4)
    ax.legend(loc="lower left", fontsize=6.4)


def main():
    if COLOR:
        plt.rcParams.update({"axes.titlesize": 8.8, "font.size": 8.2})
        figsize = (FULL, 2.55)
    else:
        figsize = (FULL, 2.4)
    fig, axs = plt.subplots(1, 3, figsize=figsize)
    panel_distortion(axs[0])
    panel_validation(axs[1])
    panel_consequence(axs[2])
    for j, ax in enumerate(axs):
        panel(ax, "abc"[j], y=0.99)
    fig.tight_layout(pad=0.5, w_pad=1.4)
    name = "fig_hero_color" if COLOR else "fig_hero"
    save(fig, os.path.join(HERE, name))


if __name__ == "__main__":
    main()
