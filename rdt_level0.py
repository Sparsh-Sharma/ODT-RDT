#!/usr/bin/env python3
r"""
Level 0 verification for the strain-coupled ODT formulation (JFM-styled output).

Closures compared for initially isotropic turbulence under constant homogeneous
mean strain A_ij:

  EXACT RDT  : integrate the spectral equations
                 dk_i/dt    = -A_ji k_j
                 dphi_ij/dt = M_in phi_nj + M_jn phi_in,  M_in = -A_in + 2 k_i k_m A_mn / k^2
               over a unit-sphere ensemble of initial wavevector directions
               (Gauss-Legendre in cos theta x uniform in phi). For initially
               isotropic turbulence the dynamics are scale-independent, so
               R_ij(t) is the spherical average of the modal covariance phi_ij.

  IP closure : dR_ij/dt = P_ij + Pi^r_ij,
               Pi^r_ij = -C2 ( P_ij - (1/3) P_kk delta_ij ),  C2 = 3/5.
               (= Sec. 3.2 baseline = LRR-IP.)

  LRR(-QI)   : Pi^r_ij = C2 k S_ij
                       + C3 k ( b_ik S_jk + b_jk S_ik - (2/3) b_mn S_mn delta_ij )
                       + C4 k ( b_ik W_jk + b_jk W_ik ),
               C2 = 0.8, C3 = 1.75, C4 = 1.31  (Launder, Reece & Rodi 1975).
               At isotropy (b=0): Pi^r = 0.8 k S = (4/5) k S  -> Crow's exact value,
               so LRR shares the exact onset slope -(8/15) k_t S_ij.

  PURE PROD. : dR_ij/dt = P_ij    (no rapid pressure-strain; shown for contrast).

Normalisation: S = (2 S_ij S_ij)^{1/2} = 1  ->  total strain e = S t = t;
R_ij(0) = (2/3) delta_ij, k_t(0) = 1.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ----------------------------------------------------------------------
# JFM / LaTeX plotting style
# ----------------------------------------------------------------------
plt.rcParams.update({
    # Computer Modern (LaTeX) look via matplotlib's bundled fonts -- no external LaTeX needed.
    # For an exact match to a LaTeX manuscript on your machine, set text.usetex=True instead.
    "text.usetex": False,
    "font.family": "serif",
    "font.serif": ["cmr10", "CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
    "axes.formatter.use_mathtext": True,
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "legend.fontsize": 8.5,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "figure.dpi": 150,
})

C2_IP = 3.0 / 5.0
C2_LRR, C3_LRR, C4_LRR = 0.8, 1.75, 1.31
I3 = np.eye(3)
EMAX, NOUT = 4.0, 200


def strain_tensor(kind):
    if kind == "plane":
        a = 0.5
        A = np.diag([a, -a, 0.0])
    elif kind == "axisymmetric":
        g = 1.0 / np.sqrt(3.0)
        A = np.diag([-0.5 * g, -0.5 * g, g])
    else:
        raise ValueError(kind)
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
    dirs = np.array(dirs); w = np.array(w); w /= w.sum()
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
        M = -A[None] + 2.0 * np.einsum("ni,nj->nij", K, AtK) / k2[:, None, None]
        Kdot = -(K @ A)
        PHIdot = np.einsum("nik,nkj->nij", M, PHI) + np.einsum("nik,njk->nij", PHI, M)
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
    sol = solve_ivp(model_rhs, [0, emax], ((2.0 / 3.0) * I3).reshape(-1),
                    t_eval=te, args=(A, rapid), rtol=1e-10, atol=1e-12)
    return te, sol.y.reshape(3, 3, -1)


def diag_fracs(R):
    kt = 0.5 * (R[0, 0] + R[1, 1] + R[2, 2])
    return [R[i, i] / (2 * kt) for i in range(3)]


# ---- compute -----------------------------------------------------------
data = {}
for kind in ["plane", "axisymmetric"]:
    A = strain_tensor(kind)
    te, Re = exact_rdt(A)
    _, Rip = integrate_model(A, "IP")
    _, Rlr = integrate_model(A, "LRR")
    _, Rpp = integrate_model(A, "none")
    data[kind] = dict(A=A, e=te,
        exact=np.array([diag_fracs(Re[:, :, t]) for t in range(NOUT)]).T,
        IP=np.array([diag_fracs(Rip[:, :, t]) for t in range(NOUT)]).T,
        LRR=np.array([diag_fracs(Rlr[:, :, t]) for t in range(NOUT)]).T,
        pp=np.array([diag_fracs(Rpp[:, :, t]) for t in range(NOUT)]).T)

print("=" * 72)
print("ONSET-SLOPE CHECK  dR_ij/de|_0  vs  -(8/15) S_ij   (k_t = 1)")
print("=" * 72)
for kind in ["plane", "axisymmetric"]:
    A = data[kind]["A"]; S = 0.5 * (A + A.T); R0 = (2.0 / 3.0) * I3
    ana = np.diag(-(8.0 / 15.0) * S)
    ip = np.diag(production(A, R0) + rapid_IP(A, R0))
    lr = np.diag(production(A, R0) + rapid_LRR(A, R0))
    print(f"\n{kind}:")
    print(f"  analytic : {np.array2string(ana, precision=5, sign='+')}")
    print(f"  IP       : {np.array2string(ip, precision=5, sign='+')}   (max err {np.max(np.abs(ip-ana)):.1e})")
    print(f"  LRR      : {np.array2string(lr, precision=5, sign='+')}   (max err {np.max(np.abs(lr-ana)):.1e})")


# ---- plotting ----------------------------------------------------------
ccol = ["#1f77b4", "#d62728", "#2ca02c"]
cnm = [r"$\overline{u_1^2}/2k_t$", r"$\overline{u_2^2}/2k_t$", r"$\overline{u_3^2}/2k_t$"]
titles = {"plane": r"Plane strain, $A=\mathrm{diag}(a,-a,0)$",
          "axisymmetric": r"Axisymmetric strain, $A=\mathrm{diag}(-\gamma/2,-\gamma/2,\gamma)$"}
mstyle = {"exact": (0, ()), "LRR": (0, (6, 2)), "IP": (0, (1, 1.5))}

def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(f"{name}.{ext}", bbox_inches="tight")
    print(f"  wrote {name}.pdf / {name}.png")

comp_sel = {"plane": [0, 1, 2], "axisymmetric": [0, 2]}
comp_lab = {"plane": {0: r"$\overline{u_1^2}/2k_t$", 1: r"$\overline{u_2^2}/2k_t$",
                      2: r"$\overline{u_3^2}/2k_t$"},
            "axisymmetric": {0: r"$\overline{u_1^2}=\overline{u_2^2}/2k_t$ (lateral)",
                             2: r"$\overline{u_3^2}/2k_t$ (axial)"}}

for kind in ["plane", "axisymmetric"]:
    d = data[kind]
    fig, ax = plt.subplots(figsize=(5.2, 3.7))
    for i in comp_sel[kind]:
        ax.plot(d["e"], d["exact"][i], color=ccol[i], ls=mstyle["exact"], lw=1.7)
        ax.plot(d["e"], d["LRR"][i],   color=ccol[i], ls=mstyle["LRR"],   lw=1.4)
        ax.plot(d["e"], d["IP"][i],    color=ccol[i], ls=mstyle["IP"],    lw=1.4)
    ax.set_xlabel(r"total strain $e = S\,t$")
    ax.set_ylabel(r"component energy fraction")
    ax.set_title(titles[kind])
    ax.set_xlim(0, EMAX)
    ax.grid(alpha=0.25, lw=0.5)
    leg_comp = ax.legend([Line2D([0], [0], color=ccol[i], lw=2) for i in comp_sel[kind]],
                         [comp_lab[kind][i] for i in comp_sel[kind]],
                         title=r"component", loc="upper left",
                         bbox_to_anchor=(1.02, 1.0), frameon=False, handlelength=1.8)
    ax.legend([Line2D([0], [0], color="k", ls=mstyle[m], lw=1.5)
               for m in ["exact", "LRR", "IP"]],
              [r"exact RDT", r"model, LRR closure", r"model, IP closure"],
              title=r"closure", loc="upper left",
              bbox_to_anchor=(1.02, 0.55), frameon=False, handlelength=2.6)
    ax.add_artist(leg_comp)
    save(fig, f"fig_components_{kind}")
    plt.close(fig)

d = data["plane"]
fig, ax = plt.subplots(figsize=(5.0, 3.7))
ax.plot(d["e"], d["exact"][1], color="k",       ls="-",                  lw=1.8, label=r"exact RDT")
ax.plot(d["e"], d["LRR"][1],   color="#1f77b4", ls=(0, (6, 2)),          lw=1.6, label=r"model, LRR closure")
ax.plot(d["e"], d["IP"][1],    color="#d62728", ls=(0, (4, 1.5, 1, 1.5)),lw=1.6, label=r"model, IP closure")
ax.plot(d["e"], d["pp"][1],    color="0.45",    ls=(0, (1, 1.5)),        lw=1.5, label=r"pure production")
ax.axhline(1.0 / 3.0, color="0.75", lw=0.7)
ax.text(0.08, 1.0 / 3.0 + 0.012, r"isotropic ($1/3$)", fontsize=8, color="0.4")
ax.set_xlabel(r"total strain $e = S\,t$")
ax.set_ylabel(r"upwash energy fraction $\overline{u_2^2}/2k_t$")
ax.set_title(r"Upwash amplification under plane strain")
ax.set_xlim(0, EMAX); ax.set_ylim(0, 1.0)
ax.grid(alpha=0.25, lw=0.5)
ax.legend(loc="lower right", frameon=False)
save(fig, "fig_upwash_amplification")
plt.close(fig)

print("\n" + "=" * 72)
print("Plane-strain upwash fraction u2^2/2kt at e = 4 :")
print(f"  exact RDT        : {data['plane']['exact'][1,-1]:.4f}")
print(f"  model (LRR)      : {data['plane']['LRR'][1,-1]:.4f}")
print(f"  model (IP)       : {data['plane']['IP'][1,-1]:.4f}")
print(f"  pure production  : {data['plane']['pp'][1,-1]:.4f}")
print("=" * 72)
