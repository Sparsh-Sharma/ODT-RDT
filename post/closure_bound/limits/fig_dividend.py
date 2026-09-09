#!/usr/bin/env python3
"""The mechanism dividend, shown in the model's own observables (for
Alan + ESM): how the size-capped relaxation clock (RCS1) improves the
strained spectra over the SC-ODT baseline, and what the frozen-von-Karman
practice misses entirely.

(a) deviation from the exact-RDT-distorted vK family at e=2, rapid
    (Sk/eps~8): binned E2 and Eperp over their best-family curves,
    standard vs clock - the fit-cost result made visible.
(b) rms log-distance from the family vs strain, both rapidities,
    standard vs clock vs the frozen-isotropic fit (the vK-practice
    surrogate, which cannot represent the strained anisotropy at all).
(c) the allocation observable at e=1, rapid: component ratio
    E2/Eperp(k2) - measured standard, measured clock, the exact-RDT
    family curve (DNS-anchored rapid reference), and the frozen
    prediction (ratio = 1 identically).

Data: spectra_gateA_{S1,S20}[_RCS1].npz + fit tables (lenoise/), 1024
realisations each, same binning/band as the fits ([3dk,300], 40 bins).
"""
import os
import re
import sys

import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
LEN = os.path.normpath(os.path.join(HERE, "..", "lenoise"))
sys.path.insert(0, LEN)
sys.path.insert(0, os.path.join(HERE, ".."))
import rdt_kernel as rk  # noqa: E402
from figstyle_jfm import BLACK, GOLD, RED, FULL, panel, plt, save  # noqa: E402,E501

NBIN, KMAX, KMIN_FAC, NKP, NPHI = 40, 300.0, 3.0, 200, 96
EOFF = {"S1": (0.4, 1.0), "S20": (0.4, 20.0)}       # (E_OFF, SMAG)


def log_bin(k, arrs, nbin):
    edges = np.geomspace(k.min(), k.max(), nbin + 1)
    idx = np.clip(np.digitize(k, edges) - 1, 0, nbin - 1)
    kb, out = [], [[] for _ in arrs]
    for b in range(nbin):
        m = idx == b
        if m.any():
            kb.append(np.exp(np.mean(np.log(k[m]))))
            for j, a in enumerate(arrs):
                out[j].append(a[m].mean())
    return np.array(kb), [np.array(o) for o in out]


def rdt_line_spectra(k2, e, ke, A0):
    kp, lnkp = rk.af._log_grid(ke, 14.0, NKP)
    phi = (np.arange(NPHI) + 0.5) * (2.0 * np.pi / NPHI)
    K1 = kp[None, :, None] * np.cos(phi)[None, None, :]
    K3 = kp[None, :, None] * np.sin(phi)[None, None, :]
    K2 = np.asarray(k2, float)[:, None, None] * np.ones_like(K1)
    P22, Ppp = rk.rdt_phi_components(K1 / ke, K2 / ke, K3 / ke, e)
    w = 2.0 * np.pi * kp * kp
    E2 = A0 * np.trapezoid(P22.mean(axis=2) * w[None, :], lnkp, axis=1)
    Ep = A0 * 0.5 * np.trapezoid(Ppp.mean(axis=2) * w[None, :], lnkp,
                                 axis=1)
    return Ep, E2


def read_fits(tag):
    fits = {}
    for line in open(os.path.join(LEN, f"rdt_family_fit_table_{tag}.txt")):
        m = re.match(r"\s*([\d.]+)\s+([\d.]+)\s+([\d.eE+-]+)\s+([\d.]+)",
                     line)
        if m:
            fits[float(m.group(1))] = (float(m.group(2)),
                                       float(m.group(3)),
                                       float(m.group(4)))
    return fits


def binned(tag, variant, j):
    npz = os.path.join(LEN, f"spectra_gateA_{tag}{variant}.npz")
    d = np.load(npz, allow_pickle=True)
    k2 = d[f"k2_e{j}"]
    Ep = 0.5 * (d[f"med_E1_e{j}"] + d[f"med_E3_e{j}"])
    E2 = d[f"med_E2_e{j}"]
    sel = (k2 >= KMIN_FAC * k2[0]) & (k2 <= KMAX)
    return log_bin(k2[sel], [Ep[sel], E2[sel]], NBIN)


def rms_from(k2b, Epb, E2b, e, ke, A0):
    Epm, E2m = rdt_line_spectra(k2b, e, ke, A0)
    res = np.concatenate([np.log(Epb / Epm), np.log(E2b / E2m)])
    return float(np.sqrt(np.mean(res ** 2)))


def fit_frozen(k2b, Epb, E2b):
    """Best ISOTROPIC (e=0) member: what a frozen-form practice can
    represent of the strained spectra; the anisotropy is pure residual."""
    ln_t = np.concatenate([np.log(Epb), np.log(E2b)])

    def resid(th):
        Ep, E2 = rdt_line_spectra(k2b, 0.0, np.exp(th[0]), np.exp(th[1]))
        return np.concatenate([np.log(np.clip(Ep, 1e-300, None)),
                               np.log(np.clip(E2, 1e-300, None))]) - ln_t

    th0 = np.array([np.log(k2b[np.argmax(E2b)]), 0.0])
    Ep0, _ = rdt_line_spectra(k2b, 0.0, np.exp(th0[0]), 1.0)
    th0[1] = np.log(np.median(Epb) / max(np.median(Ep0), 1e-300))
    sol = least_squares(resid, th0, method="trf", max_nfev=60)
    return np.exp(sol.x[0]), np.exp(sol.x[1])


def main():
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(FULL, 2.35))

    # ---- (a) deviation from the family, S20 e=2 ------------------------
    j, e = 4, 2.0
    for variant, col, lab in ((" ", BLACK, "SC-ODT"),
                              ("_RCS1", RED, "clock-relaxed ODT")):
        v = "" if variant == " " else variant
        k2b, (Epb, E2b) = binned("S20", v, j)
        ke, A0, _ = read_fits("S20" + v)[e]
        Epm, E2m = rdt_line_spectra(k2b, e, ke, A0)
        ke0 = read_fits("S20" + v)[0.0][0]
        a.semilogx(k2b / ke0, E2b / E2m, "-", color=col, lw=1.0,
                   label=lab)
        a.semilogx(k2b / ke0, Epb / Epm, ls=(0, (2, 1.5)), color=col,
                   lw=0.8, alpha=0.7)
    a.axhline(1.0, color="0.6", lw=0.6, ls=":")
    a.set_xlabel(r"$\kappa_2/\kappa_e(0)$")
    a.set_ylabel(r"$E_i/E_i^{\rm RDT\,fit}$")
    a.set_xticks([6, 10, 20, 40])
    a.set_xticklabels(["6", "10", "20", "40"])
    a.minorticks_off()
    a.legend(fontsize=6, loc="upper right")

    # ---- (b) rms distance vs strain, both rapidities -------------------
    for tag, mk in (("S1", "o"), ("S20", "s")):
        es = [0.5, 1.0, 1.5, 2.0]
        for v, col in (("", BLACK), ("_RCS1", RED)):
            fits = read_fits(tag + v)
            r = []
            for e in es:
                jj = int(round(e / 0.5))
                k2b, (Epb, E2b) = binned(tag, v, jj)
                ke, A0, _ = fits[e]
                r.append(100 * rms_from(k2b, Epb, E2b, e, ke, A0))
            b.plot(es, r, marker=mk, ms=3, color=col, lw=0.9,
                   ls="-" if tag == "S20" else (0, (4, 2)))
        # frozen-isotropic surrogate on the STANDARD ensembles
        r = []
        for e in es:
            jj = int(round(e / 0.5))
            k2b, (Epb, E2b) = binned(tag, "", jj)
            ke, A0 = fit_frozen(k2b, Epb, E2b)
            r.append(100 * rms_from(k2b, Epb, E2b, 0.0, ke, A0))
        b.plot(es, r, marker=mk, ms=3, color=GOLD, lw=0.9,
               ls="-" if tag == "S20" else (0, (4, 2)))
    b.set_xlabel(r"total strain $e$")
    b.set_ylabel(r"rms distance from RDT family [\%]")
    from matplotlib.lines import Line2D
    b.legend(handles=[
        Line2D([], [], color=GOLD, lw=1, label="frozen isotropic"),
        Line2D([], [], color=BLACK, lw=1, label="SC-ODT"),
        Line2D([], [], color=RED, lw=1, label="clock-relaxed ODT"),
        Line2D([], [], color="0.4", marker="s", ls="-", ms=3, lw=0.8,
               label=r"$Sk_t/\varepsilon\approx8$"),
        Line2D([], [], color="0.4", marker="o", ls=(0, (4, 2)), ms=3,
               lw=0.8, label=r"$\approx0.4$")],
        fontsize=5.5, loc="upper left", ncol=2, columnspacing=0.8)

    # ---- (c) allocation observable, S20 e=1 ----------------------------
    j, e = 2, 1.0
    ke_r, A0_r, _ = read_fits("S20_RCS1")[e]
    ke0 = read_fits("S20_RCS1")[0.0][0]
    for v, col, lab in (("", BLACK, "SC-ODT"),
                        ("_RCS1", RED, "clock-relaxed ODT")):
        k2b, (Epb, E2b) = binned("S20", v, j)
        c.semilogx(k2b / ke0, E2b / Epb, "-", color=col, lw=1.0,
                   label=lab)
    kk = np.geomspace(KMIN_FAC * 0.4, KMAX, 60)
    Epm, E2m = rdt_line_spectra(kk, e, ke_r, A0_r)
    c.semilogx(kk / ke0, E2m / Epm, color="#8e44ad", lw=1.2,
               ls=(0, (5, 2)), label="exact RDT")
    c.axhline(1.0, color=GOLD, lw=1.0, ls=":", label="frozen (any form)")
    c.set_xlabel(r"$\kappa_2/\kappa_e(0)$")
    c.set_ylabel(r"$E_2/E_\perp$")
    c.set_ylim(0.8, 2.6)
    c.text(0.05, 0.62, "RDT ratio $\\to8$ at the\nenergy scales "
           "(off-scale)", transform=c.transAxes, fontsize=5.5,
           color="#8e44ad")
    c.legend(fontsize=6, loc="upper right")

    for ax, s in ((a, "a"), (b, "b"), (c, "c")):
        panel(ax, s)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_dividend"))


if __name__ == "__main__":
    main()
