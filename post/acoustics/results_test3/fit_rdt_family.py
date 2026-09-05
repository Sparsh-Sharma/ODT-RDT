#!/usr/bin/env python3
"""Gate A v2: fit the RDT-DISTORTED von Karman family to the gateA_S1
median line spectra.

Instead of the axisymmetric stretch ansatz (whose single-argument A-shape
exact RDT violates — the h(mu) swap alone left the fits degenerate), the
strained model is the EXACT plane-strain Cauchy-RDT map applied to an
isotropic vK spectrum:  two free parameters per strain, (A0, ke), the
pre-distortion amplitude and energy wavenumber; e is known.

    E_2(k2; e)    = int 2 pi kp <Phi_22^RDT>_ring dkp
    E_perp(k2; e) = int 2 pi kp <(Phi_11+Phi_33)/2^RDT>_ring dkp

both computed numerically from rdt_kernel.rdt_phi_components on a log-kp
grid.  Log-space least squares on (Eperp, E2), same band/binning as the
axisym fits.

    E_OFF=0.4 KMAX=600 python3 fit_rdt_family.py spectra_gateA_S1.npz
"""
import os
import sys

import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import odt_io                                          # noqa: E402
import rdt_kernel as rk                                # noqa: E402

NPZ = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    HERE, "spectra_gateA_S1.npz")
NBIN = 40
KMAX = float(os.environ.get("KMAX", 600.0))
KMIN_FAC = float(os.environ.get("KMIN_FAC", 3.0))
E_OFF = float(os.environ.get("E_OFF", 0.4))
SMAG = float(os.environ.get("SMAG", 1.0))    # e = SMAG * (dump time - E_OFF)
NKP = 200
NPHI = 96


def rdt_line_spectra(k2_grid, e, ke, A0):
    """Model (Eperp, E2) on k2_grid for the RDT-distorted vK(ke, A0)."""
    k2 = np.asarray(k2_grid, float)
    kp, lnkp = rk.af._log_grid(ke, 14.0, NKP)
    phi = (np.arange(NPHI) + 0.5) * (2.0 * np.pi / NPHI)
    K1 = kp[None, :, None] * np.cos(phi)[None, None, :]
    K3 = kp[None, :, None] * np.sin(phi)[None, None, :]
    K2 = k2[:, None, None] * np.ones_like(K1)
    # vK scale enters through E0(k0/ke): evaluate with scaled wavevectors
    P22, Ppp = rk.rdt_phi_components(K1 / ke, K2 / ke, K3 / ke, e)
    # E0 argument scaling: Phi built from E0(|k0|) with ke=1; k -> k/ke and
    # an overall ke^-3 Jacobian are absorbed into A0 (fit amplitude).
    r22 = P22.mean(axis=2)
    rpp = Ppp.mean(axis=2)
    w = 2.0 * np.pi * kp * kp                       # int f dkp = int f kp dln
    E2 = A0 * np.trapezoid(r22 * w[None, :], lnkp, axis=1)
    Ep = A0 * 0.5 * np.trapezoid(rpp * w[None, :], lnkp, axis=1)
    return Ep, E2


def fit_one(k2b, Epb, E2b, e):
    ln_t = np.concatenate([np.log(Epb), np.log(E2b)])

    def resid(theta):
        lnke, lnA0 = theta
        Ep, E2 = rdt_line_spectra(k2b, e, np.exp(lnke), np.exp(lnA0))
        return np.concatenate([np.log(np.clip(Ep, 1e-300, None)),
                               np.log(np.clip(E2, 1e-300, None))]) - ln_t

    kpk = k2b[np.argmax(E2b)]
    th0 = np.array([np.log(max(kpk, 1.0)), 0.0])
    # amplitude: one cheap evaluation to center A0
    Ep0, _ = rdt_line_spectra(k2b, e, np.exp(th0[0]), 1.0)
    th0[1] = np.log(np.median(Epb) / max(np.median(Ep0), 1e-300))
    sol = least_squares(resid, th0, method="trf", xtol=1e-10, ftol=1e-10,
                        max_nfev=80)
    return np.exp(sol.x[0]), np.exp(sol.x[1]), sol.cost, sol


def main():
    d = np.load(NPZ)
    strains = SMAG * (d["strains"] - E_OFF)
    lines = [f"RDT-distorted vK fits, median spectra, {int(d['nrlz'])} rlz "
             f"({os.path.basename(NPZ)}); band [3dk,{KMAX:.0f}], {NBIN} bins",
             "", f"{'e':>5} {'ke':>9} {'A0':>11} {'cost':>9}"]
    rows = []
    for j, e in enumerate(strains):
        k2 = d[f"k2_e{j}"]
        E1, E2, E3 = d[f"med_E1_e{j}"], d[f"med_E2_e{j}"], d[f"med_E3_e{j}"]
        Eperp = 0.5 * (E1 + E3)
        sel = (k2 >= KMIN_FAC * k2[0]) & (k2 <= KMAX)
        k2b, (Epb, E2b) = odt_io.log_bin(k2[sel], [Eperp[sel], E2[sel]], NBIN)
        ke, A0, cost, _ = fit_one(k2b, Epb, E2b, e)
        rows.append((e, ke, A0, cost, k2b, Epb, E2b))
        lines.append(f"{e:5.1f} {ke:9.3f} {A0:11.4e} {cost:9.3f}")
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(HERE, "rdt_family_fit_table.txt"), "w").write(txt + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.4))
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(rows)))
    for (e, ke, A0, cost, k2b, Epb, E2b), col in zip(rows, cols):
        Epm, E2m = rdt_line_spectra(k2b, e, ke, A0)
        axs[0].loglog(k2b, E2b, "o", ms=2.5, color=col)
        axs[0].loglog(k2b, E2m, "-", lw=1.1, color=col,
                      label=f"e={e:.1f}: ke={ke:.1f}, cost={cost:.2f}")
        axs[1].loglog(k2b, Epb, "o", ms=2.5, color=col)
        axs[1].loglog(k2b, Epm, "-", lw=1.1, color=col)
    axs[0].set_title("$E_2$: data + RDT-vK fit")
    axs[1].set_title(r"$E_\perp$: data + RDT-vK fit")
    for ax in axs:
        ax.set_xlabel("$k_2$")
    axs[0].legend(fontsize=7)
    fig.suptitle("Gate A v2: exact-RDT-distorted vK (2 params/strain)",
                 fontsize=11)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"fig_rdt_family.{ext}"), dpi=200)
    print("saved fig_rdt_family.png/.pdf")


if __name__ == "__main__":
    main()
