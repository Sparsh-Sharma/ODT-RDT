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

def load(tag):
    return np.load(os.path.join(HERE, f"{tag}_ens.npz"))


def fig_bvc(on, off):
    # monochrome: SC-ODT (strain on) = solid + grey filled band;
    # standard ODT (strain off) = dashed + hatched open band.
    ON = dict(color="k", ls="-", lw=1.2)
    OFF = dict(color="k", ls=(0, (5, 2)), lw=1.0)
    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.3))
    k, E, gs = on["k_last"], on["Etot_last"], on["Etot_last_gsd"]
    a.fill_between(k, E / gs, E * gs, color="0.8", alpha=0.6, lw=0)
    a.loglog(k, E, label="SC-ODT (strain on)", **ON)
    k, E, gs = off["k_last"], off["Etot_last"], off["Etot_last_gsd"]
    a.fill_between(k, E / gs, E * gs, facecolor="none", edgecolor="0.55",
                   hatch="////", lw=0.0)
    a.loglog(k, E, label="standard ODT (strain off)", **OFF)
    a.set_xlabel(r"wavenumber $\kappa_2$")
    a.set_ylabel(r"$E(\kappa_2)$")
    a.set_xlim(6, 2e4)
    ymax = on["Etot_last"].max() * 3
    a.set_ylim(ymax * 1e-9, ymax)
    a.legend(loc="lower left")
    ee = np.linspace(0, 4, 100)
    b.plot(ee, np.exp(0.5 * ee), color="0.45", ls=(0, (1, 1.2)), lw=1.0,
           label="rigid translation")
    for d, st, lab in ((on, ON, "SC-ODT (strain on)"),
                       (off, OFF, "standard ODT (strain off)")):
        r = d["cent_gmean"] / d["cent_gmean"][0]
        gsd = d["cent_gsd"]
        b.plot(d["t"], r, label=lab, **st)
        if lab.startswith("SC"):
            b.fill_between(d["t"], r / gsd, r * gsd, color="0.8",
                           alpha=0.6, lw=0)
        else:
            b.fill_between(d["t"], r / gsd, r * gsd, facecolor="none",
                           edgecolor="0.55", hatch="////", lw=0.0)
    b.axhline(1.0, color="0.7", lw=0.5, ls=":")
    b.set_xlabel(r"total strain $e$")
    b.set_ylabel(r"$\bar\kappa(e)/\bar\kappa(0)$")
    b.set_yscale("log")
    b.legend(loc="upper left")
    panel(a, "a", x=0.88)
    panel(b, "b", x=0.88, y=0.15)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_spec_bvc"))
    print("bvc: strain-on centroid ratio at e=3.9 (vs t=0):",
          f"{on['cent_gmean'][-1] / on['cent_gmean'][0]:.2f}")
    print("bvc: strain-off centroid ratio at e=3.9 (vs t=0):",
          f"{off['cent_gmean'][-1] / off['cent_gmean'][0]:.2f}")


def fig_op(op):
    # NB: until the operating-point caro rerun (nu=3e-6, N=200,
    # input/homogeneousStrainOP) lands as a real OP_ens.npz, the
    # manuscript figure fig_spec_op is built by fig_spec_op_archive.py
    # from the archived vectors -- NOT here.  This path is guarded so a
    # degenerate (single-timestep) ensemble cannot silently overwrite it.
    if np.ndim(op["t"]) == 0 or len(op["t"]) < 2:
        raise ValueError(
            "OP ensemble has < 2 timesteps (stub/degenerate); refusing to "
            "regenerate fig_spec_op. Use fig_spec_op_archive.py, or supply "
            "a real OP_ens.npz from the caro rerun.")
    # monochrome: components -> line style (E1 solid, E2 dashed, E3
    # dotted), matching fig_spec_kinematic.
    LS = ("-", (0, (6, 2)), (0, (1, 1.4)))
    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.3))
    for i, lab in enumerate(("$E_1$", "$E_2$", "$E_3$")):
        a.loglog(op["k_last"], op[f"E{i + 1}_last"], color="k",
                 ls=LS[i], lw=1.0, label=lab)
    kref = np.array([200.0, 4000.0])
    a.plot(kref, 3e-3 * op["Etot_last"].max()
           * (kref / kref[0]) ** (-5 / 3), color="0.45", ls=(0, (1, 1.2)),
           lw=0.9)
    a.text(1000, 4e-3 * op["Etot_last"].max()
           * (1000 / kref[0]) ** (-5 / 3), r"$\kappa^{-5/3}$",
           fontsize=7, color="0.3")
    a.set_xlabel(r"wavenumber $\kappa_2$")
    a.set_ylabel(r"$E_i(\kappa_2)$")
    ymax = op["Etot_last"].max() * 3
    a.set_xlim(20, 1e5)
    a.set_ylim(ymax * 1e-10, ymax)
    a.legend(loc="lower left")
    ee = np.linspace(0, 4, 100)
    b.plot(ee, np.exp(0.5 * ee), color="0.45", ls=(0, (1, 1.2)), lw=1.0,
           label="rigid translation")
    r = op["cent_gmean"] / op["cent_gmean"][0]
    b.plot(op["t"], r, color="k", ls="-", lw=1.2, label="SC-ODT centroid")
    b.fill_between(op["t"], r / op["cent_gsd"], r * op["cent_gsd"],
                   color="0.8", alpha=0.6, lw=0)
    rp = op["peak_gmean"] / op["peak_gmean"][0]
    b.plot(op["t"], rp, color="k", lw=1.0, ls=(0, (5, 2)),
           label="spectral peak")
    b.set_xlabel(r"total strain $e$")
    b.set_ylabel(r"migration from $e=0$")
    b.set_yscale("log")
    b.legend(loc="upper left")
    panel(a, "a", x=0.88)
    panel(b, "b", x=0.88, y=0.15)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_spec_op"))
    print("op: centroid ratio at e=3.9:", f"{r[-1]:.2f}",
          " peak ratio:", f"{rp[-1]:.2f}",
          " rigid:", f"{np.exp(0.5 * op['t'][-1]):.2f}")
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


def fig_aniso(on, off):
    """Strain-induced migration per component: centroid of the strained
    ensemble over the matched no-strain baseline (the on/off logic of
    sec:spec-bvc applied per component), against the rigid law."""
    fig, ax = plt.subplots(figsize=(0.62 * FULL, 2.3))
    t = on["t"]
    ee = np.linspace(0, 2.2, 50)
    ax.plot(ee, np.exp(0.5 * ee), "k--", lw=0.9,
            label="rigid translation")
    for key, col, lab in (("c1", COL["phi11"],
                           r"streamwise $\bar\kappa_1$"),
                          ("c2", COL["phi22"],
                           r"upwash $\bar\kappa_2$")):
        r = on[key + "_gmean"] / off[key + "_gmean"]
        gs = np.sqrt(on[key + "_gsd"] * off[key + "_gsd"])
        ax.plot(t, r, color=col, lw=1.1, label=lab)
        ax.fill_between(t, r / gs, r * gs, color=col, alpha=0.2, lw=0)
    ax.axhline(1.0, color="0.7", lw=0.5, ls=":")
    ax.set_xlim(0, 2.2)
    ax.set_ylim(0.9, 3.2)
    ax.set_xlabel(r"total strain $e$")
    ax.set_ylabel(r"$\bar\kappa_i^{\rm on}/\bar\kappa_i^{\rm off}$")
    ax.legend(loc="upper left")
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_cmk_aniso"))
    for key in ("c1", "c2"):
        r = on[key + "_gmean"] / off[key + "_gmean"]
        print(f"aniso: {key} on/off at e=2.2:",
              f"{np.interp(2.2, t, r):.2f}")


if __name__ == "__main__":
    on = load("B5")
    off = load("B5off")
    fig_bvc(on, off)
    fig_aniso(on, off)
    try:
        op = load("OP")
    except FileNotFoundError:
        print("OP ensemble not ready (caro rerun pending); "
              "fig_spec_op is built by fig_spec_op_archive.py")
    else:
        fig_op(op)
