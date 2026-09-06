#!/usr/bin/env python3
"""Manuscript figures for sec:verification (main_v3.tex): moment-level
verification of the strain-coupled formulation against exact RDT.

Computation (from the original Level-0 script rdt_level0.py, verbatim):
exact RDT by unit-sphere wavevector-ensemble integration; the model
moment equation with the IP (C2=3/5) and LRR-QI closures; pure
production for contrast.  S=1 so e=t; R(0)=(2/3)I.

Outputs (JFM canon, figstyle_jfm):
  fig_rdt_components.pdf : (a) plane strain, (b) axisymmetric strain —
                           component fractions, exact/LRR/IP.
  fig_upwash_amplification.pdf : plane-strain upwash fraction, incl.
                           production-only.
"""
import os
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from figstyle_jfm import COL, FULL, GOLD, panel, plt, save  # noqa: E402

C2_IP = 3.0 / 5.0
C2_LRR, C3_LRR, C4_LRR = 0.8, 1.75, 1.31
I3 = np.eye(3)
EMAX, NOUT = 4.0, 200


def strain_tensor(kind):
    if kind == "plane":
        A = np.diag([0.5, -0.5, 0.0])
    else:
        g = 1.0 / np.sqrt(3.0)
        A = np.diag([-0.5 * g, -0.5 * g, g])
    S = 0.5 * (A + A.T)
    assert abs(np.sqrt(2.0 * np.sum(S * S)) - 1.0) < 1e-12
    return A


def sphere_quadrature(nmu=24, nphi=48):
    mu, wmu = leggauss(nmu)
    phi = (np.arange(nphi) + 0.5) * 2.0 * np.pi / nphi
    dirs, w = [], []
    for m, wm in zip(mu, wmu):
        st = np.sqrt(max(0.0, 1.0 - m * m))
        for p in phi:
            dirs.append([st * np.cos(p), st * np.sin(p), m])
            w.append(0.5 * wm / nphi)
    dirs = np.array(dirs)
    w = np.array(w)
    w /= w.sum()
    return dirs, w


def exact_rdt(A, emax=EMAX, nout=NOUT, nmu=24, nphi=48):
    dirs, w = sphere_quadrature(nmu, nphi)
    N = dirs.shape[0]
    PHI0 = I3[None] - np.einsum("ni,nj->nij", dirs, dirs)
    y0 = np.concatenate([dirs.reshape(-1), PHI0.reshape(-1)])

    def rhs(t, y):
        K = y[:3 * N].reshape(N, 3)
        PHI = y[3 * N:].reshape(N, 3, 3)
        AtK = K @ A
        k2 = np.einsum("ni,ni->n", K, K)
        M = -A[None] + 2.0 * np.einsum("ni,nj->nij", K, AtK) \
            / k2[:, None, None]
        Kdot = -(K @ A)
        PHIdot = np.einsum("nik,nkj->nij", M, PHI) \
            + np.einsum("nik,njk->nij", PHI, M)
        return np.concatenate([Kdot.reshape(-1), PHIdot.reshape(-1)])

    te = np.linspace(0, emax, nout)
    sol = solve_ivp(rhs, [0, emax], y0, t_eval=te, rtol=1e-9, atol=1e-11)
    PHI = sol.y[3 * N:].reshape(N, 3, 3, -1)
    return te, np.einsum("n,nijt->ijt", w, PHI)


def production(A, R):
    return -(A @ R + R @ A.T)


def rapid_IP(A, R):
    P = production(A, R)
    return -C2_IP * (P - np.trace(P) / 3.0 * I3)


def rapid_LRR(A, R):
    kt = 0.5 * np.trace(R)
    b = R / (2 * kt) - I3 / 3.0
    S = 0.5 * (A + A.T)
    W = 0.5 * (A - A.T)
    t2 = b @ S + S @ b - (2.0 / 3.0) * np.trace(b @ S) * I3
    t3 = W @ b - b @ W
    return C2_LRR * kt * S + C3_LRR * kt * t2 + C4_LRR * kt * t3


def model_rhs(t, y, A, rapid):
    R = y.reshape(3, 3)
    P = production(A, R)
    Pir = {"IP": rapid_IP, "LRR": rapid_LRR,
           "none": lambda A, R: np.zeros((3, 3))}[rapid](A, R)
    return (P + Pir).reshape(-1)


def integrate_model(A, rapid, emax=EMAX, nout=NOUT):
    te = np.linspace(0, emax, nout)
    sol = solve_ivp(model_rhs, [0, emax],
                    ((2.0 / 3.0) * I3).reshape(-1),
                    t_eval=te, args=(A, rapid), rtol=1e-10, atol=1e-12)
    return te, sol.y.reshape(3, 3, -1)


def diag_fracs(R):
    kt = 0.5 * (R[0, 0] + R[1, 1] + R[2, 2])
    return [R[i, i] / (2 * kt) for i in range(3)]


def main():
    data = {}
    for kind in ("plane", "axisymmetric"):
        A = strain_tensor(kind)
        te, Re = exact_rdt(A)
        _, Rip = integrate_model(A, "IP")
        _, Rlr = integrate_model(A, "LRR")
        _, Rpp = integrate_model(A, "none")
        data[kind] = dict(
            A=A, e=te,
            exact=np.array([diag_fracs(Re[:, :, t])
                            for t in range(NOUT)]).T,
            IP=np.array([diag_fracs(Rip[:, :, t])
                         for t in range(NOUT)]).T,
            LRR=np.array([diag_fracs(Rlr[:, :, t])
                          for t in range(NOUT)]).T,
            pp=np.array([diag_fracs(Rpp[:, :, t])
                         for t in range(NOUT)]).T)

    for kind in ("plane", "axisymmetric"):
        A = data[kind]["A"]
        S = 0.5 * (A + A.T)
        R0 = (2.0 / 3.0) * I3
        ana = np.diag(-(8.0 / 15.0) * S)
        for nm, fn in (("IP", rapid_IP), ("LRR", rapid_LRR)):
            err = np.max(np.abs(np.diag(production(A, R0) + fn(A, R0))
                                - ana))
            print(f"{kind:13s} onset {nm}: max err {err:.1e}")
    d = data["plane"]
    print("plane upwash at e=4: exact %.3f LRR %.3f IP %.3f prod %.3f"
          % tuple(d[m][1, -1] for m in ("exact", "LRR", "IP", "pp")))
    da = data["axisymmetric"]
    print("axisym at e=4: lateral exact %.3f LRR %.3f IP %.3f | "
          "axial exact %.3f LRR %.3f IP %.3f"
          % (da["exact"][0, -1], da["LRR"][0, -1], da["IP"][0, -1],
             da["exact"][2, -1], da["LRR"][2, -1], da["IP"][2, -1]))

    CC = [COL["phi11"], COL["phi22"], COL["phi33"]]
    LS = {"exact": "-", "LRR": (0, (6, 2)), "IP": (0, (1, 1.5))}
    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.3), sharex=True)
    for i in range(3):
        for m in ("exact", "LRR", "IP"):
            a.plot(data["plane"]["e"], data["plane"][m][i],
                   color=CC[i], ls=LS[m],
                   lw=1.3 if m == "exact" else 1.0)
    for i in (0, 2):
        for m in ("exact", "LRR", "IP"):
            b.plot(data["axisymmetric"]["e"], data["axisymmetric"][m][i],
                   color=CC[i], ls=LS[m],
                   lw=1.3 if m == "exact" else 1.0)
    for j, ax in enumerate((a, b)):
        ax.set_xlabel(r"total strain $e = S\,t$")
        ax.set_xlim(0, EMAX)
        panel(ax, "ab"[j])
    a.set_ylabel(r"$\overline{u_i^2}/2k_t$")
    from matplotlib.lines import Line2D
    hand = ([Line2D([0], [0], color=CC[i], lw=1.5) for i in range(3)]
            + [Line2D([0], [0], color="k", ls=LS[m], lw=1.2)
               for m in ("exact", "LRR", "IP")])
    labs = [r"$\overline{u_1^2}$", r"$\overline{u_2^2}$",
            r"$\overline{u_3^2}$",
            "exact RDT", "LRR closure", "IP closure"]
    fig.legend(hand, labs, ncol=6, loc="lower center",
               bbox_to_anchor=(0.5, -0.13), columnspacing=1.2)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_rdt_components"))

    fig, ax = plt.subplots(figsize=(0.62 * FULL, 2.3))
    ax.plot(d["e"], d["exact"][1], color=COL["rdt"], lw=1.4,
            label="exact RDT")
    ax.plot(d["e"], d["LRR"][1], color=COL["odt"], ls=(0, (6, 2)),
            label="model, LRR closure")
    ax.plot(d["e"], d["IP"][1], color=GOLD, ls=(0, (4, 1.5, 1, 1.5)),
            label="model, IP closure")
    ax.plot(d["e"], d["pp"][1], color="0.45", ls=(0, (1, 1.5)),
            label="production alone")
    ax.axhline(1.0 / 3.0, color="0.75", lw=0.6)
    ax.text(0.1, 1.0 / 3.0 + 0.02, "isotropic (1/3)", fontsize=7,
            color="0.4")
    ax.set_xlabel(r"total strain $e = S\,t$")
    ax.set_ylabel(r"$\overline{u_2^2}/2k_t$")
    ax.set_xlim(0, EMAX)
    ax.set_ylim(0, 1.0)
    ax.legend(loc="upper left", fontsize=6.5)
    fig.tight_layout(pad=0.4)
    save(fig, os.path.join(HERE, "fig_upwash_amplification"))


if __name__ == "__main__":
    main()
