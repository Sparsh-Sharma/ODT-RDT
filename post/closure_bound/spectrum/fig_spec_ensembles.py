#!/usr/bin/env python3
"""Manuscript figures for sec:spec-bvc / sec:spec-op / sec:cmk-full
(main_v3.tex), from the caro ensemble campaigns of 2026-09-07
(cases homogeneousStrainB5 / B5off / OP; postprocessed by
ens_spectra.py into *_ens.npz):

  fig_spec_bvc.pdf  : (a) final-strain total spectra, strain on vs off;
                      (b) centroid migration vs the rigid law.
  fig_spec_op.pdf   : (a) final component spectra at the operating
                      point; (b) total-centroid migration vs rigid law.
  fig_cmk_aniso.pdf : component-centroid divergence (E1 vs E2) of the
                      nu=1e-5 strained ensemble, normalised to the
                      relaxed state after the initialisation transient.

Also prints the numbers quoted in the text (centroid plateaus, chi from
the energy budget, peak migration).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from figstyle_jfm import COL, FULL, panel, plt, save  # noqa: E402

E0 = 0.5    # relaxed-state reference strain for normalisations


def load(tag):
    return np.load(os.path.join(HERE, f"{tag}_ens.npz"))


def jref(d, e):
    return int(np.argmin(np.abs(d["t"] - e)))


def fig_bvc(on, off):
    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.3))
    for d, col, lab in ((on, COL["odt"], "strain on"),
                        (off, "0.35", "strain off")):
        k, E, gs = d["k_last"], d["Etot_last"], d["Etot_last_gsd"]
        a.loglog(k, E, color=col, lw=1.1, label=lab)
        a.fill_between(k, E / gs, E * gs, color=col, alpha=0.2, lw=0)
    a.set_xlabel(r"wavenumber $\kappa_2$")
    a.set_ylabel(r"$E(\kappa_2)$")
    a.set_xlim(6, 2e4)
    ymax = on["Etot_last"].max() * 3
    a.set_ylim(ymax * 1e-9, ymax)
    a.legend(loc="lower left")
    ee = np.linspace(0, 4, 100)
    b.plot(ee, np.exp(0.5 * ee), "k--", lw=0.9,
           label="rigid translation")
    for d, col, lab in ((on, COL["odt"], "strain on"),
                        (off, "0.35", "strain off")):
        j0 = jref(d, E0)
        r = d["cent_gmean"] / d["cent_gmean"][j0]
        gs = d["cent_gsd"]
        b.plot(d["t"], r, color=col, lw=1.1, label=lab)
        b.fill_between(d["t"], r / gs, r * gs, color=col, alpha=0.2,
                       lw=0)
    b.axhline(1.0, color="0.7", lw=0.5, ls=":")
    b.set_xlabel(r"total strain $e$")
    b.set_ylabel(r"$\bar\kappa(e)/\bar\kappa(e_0)$")
    b.set_yscale("log")
    b.legend(loc="upper left")
    panel(a, "a", x=0.88)
    panel(b, "b", x=0.88, y=0.15)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_spec_bvc"))
    j0 = jref(on, E0)
    print("bvc: strain-on centroid ratio at e=4:",
          f"{on['cent_gmean'][-1] / on['cent_gmean'][j0]:.2f}")
    j0 = jref(off, E0)
    print("bvc: strain-off centroid ratio at e=4:",
          f"{off['cent_gmean'][-1] / off['cent_gmean'][j0]:.2f}")


def fig_op(op):
    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.3))
    for i, (col, lab) in enumerate(
            ((COL["phi11"], "$E_1$"), (COL["phi22"], "$E_2$"),
             (COL["phi33"], "$E_3$"))):
        a.loglog(op["k_last"], op[f"E{i + 1}_last"], color=col, lw=0.9,
                 label=lab)
    kref = np.array([200.0, 4000.0])
    a.plot(kref, 3e-3 * op["Etot_last"].max()
           * (kref / kref[0]) ** (-5 / 3), "k:", lw=0.7)
    a.text(1000, 4e-3 * op["Etot_last"].max()
           * (1000 / kref[0]) ** (-5 / 3), r"$\kappa^{-5/3}$",
           fontsize=7)
    a.set_xlabel(r"wavenumber $\kappa_2$")
    a.set_ylabel(r"$E_i(\kappa_2)$")
    ymax = op["Etot_last"].max() * 3
    a.set_xlim(20, 1e5)
    a.set_ylim(ymax * 1e-10, ymax)
    a.legend(loc="lower left")
    ee = np.linspace(0, 4, 100)
    b.plot(ee, np.exp(0.5 * ee), "k--", lw=0.9,
           label="rigid translation")
    j0 = jref(op, E0)
    r = op["cent_gmean"] / op["cent_gmean"][j0]
    b.plot(op["t"], r, color=COL["odt"], lw=1.1, label="ODT centroid")
    b.fill_between(op["t"], r / op["cent_gsd"], r * op["cent_gsd"],
                   color=COL["odt"], alpha=0.2, lw=0)
    rp = op["peak_gmean"] / op["peak_gmean"][j0]
    b.plot(op["t"], rp, color=COL["phi11"], lw=0.9, ls=(0, (4, 2)),
           label="spectral peak")
    b.set_xlabel(r"total strain $e$")
    b.set_ylabel(r"migration from $e_0=0.5$")
    b.set_yscale("log")
    b.legend(loc="upper left")
    panel(a, "a", x=0.88)
    panel(b, "b", x=0.88, y=0.15)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_spec_op"))
    print("op: centroid ratio at e=4:", f"{r[-1]:.2f}",
          " peak ratio:", f"{rp[-1]:.2f}",
          " rigid:", f"{np.exp(0.5 * (op['t'][-1] - E0)):.2f}")
    # chi from the energy budget:  P = 0.5*(R22 - R11) for A=diag(.5,-.5,0)
    kt, t = op["kt_mean"], op["t"]
    P = 0.5 * (op["R22_mean"] - op["R11_mean"])
    dkt = np.gradient(kt, t)
    eps = P - dkt
    chi = kt / np.maximum(eps, 1e-30)
    m = (t > 0.5) & (eps > 0)
    print("op: chi (budget) median over e>0.5:",
          f"{np.median(chi[m]):.2f}",
          " IQR", f"{np.percentile(chi[m], 25):.2f}",
          f"{np.percentile(chi[m], 75):.2f}")
    nu = 3.0e-6
    chig = kt / (nu * op["eps_g_mean"])
    print("op: chi (gradient) median:", f"{np.median(chig[m]):.2f}")


def fig_aniso(on):
    fig, ax = plt.subplots(figsize=(0.62 * FULL, 2.3))
    j0 = jref(on, E0)
    for key, col, lab in (("c1", COL["phi11"],
                           r"streamwise $\bar\kappa_1$"),
                          ("c2", COL["phi22"],
                           r"upwash $\bar\kappa_2$")):
        r = on[key + "_gmean"] / on[key + "_gmean"][j0]
        ax.plot(on["t"], r, color=col, lw=1.1, label=lab)
        ax.fill_between(on["t"], r / on[key + "_gsd"],
                        r * on[key + "_gsd"], color=col, alpha=0.2,
                        lw=0)
    ax.axvspan(0, E0, color="0.92", zorder=0)
    ax.axhline(1.0, color="0.7", lw=0.5, ls=":")
    ax.set_xlabel(r"total strain $e$")
    ax.set_ylabel(r"$\bar\kappa_i(e)/\bar\kappa_i(e_0)$")
    ax.legend(loc="upper left")
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_cmk_aniso"))
    print("aniso: c1, c2 ratio at e=2.2:",
          f"{np.interp(2.2, on['t'], on['c1_gmean'] / on['c1_gmean'][j0]):.2f}",
          f"{np.interp(2.2, on['t'], on['c2_gmean'] / on['c2_gmean'][j0]):.2f}")


if __name__ == "__main__":
    on = load("B5")
    off = load("B5off")
    fig_bvc(on, off)
    fig_aniso(on)
    try:
        op = load("OP")
    except FileNotFoundError:
        print("OP ensemble not ready; skipping fig_spec_op")
    else:
        fig_op(op)
