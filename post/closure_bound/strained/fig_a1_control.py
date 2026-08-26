"""Payoff figure: the A1 closure error is a strain-geometry effect.

Two panels vs accumulated strain e at S k/eps = 16 (rapid, effect clearest):
 (left)  azimuthal m=2 residue of Phi_22 -- the direct A1-violation measure;
 (right) transverse line-spectrum splitting phi_11/phi_33 - 1.
Plane strain (A1 broken) vs axisymmetric contraction (A1 exact). A third
line shows b_22(e), identical in both geometries -- the acoustically
relevant upwash amplification is robust; only the azimuthal structure that
A1 assumes away differs.
"""

import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ECHECKS = ("0", "0.25", "0.5", "0.75", "1")
EVAL = [float(e) for e in ECHECKS]
KMAX_FRAC = 0.85
RATIO = "16"

C_PLANE, C_AXI, C_B22 = "#eb6834", "#2a78d6", "#1baf7a"


def group(subdir, prefix, e):
    pat = os.path.join(HERE, subdir, f"chk_{prefix}r{RATIO}_s*_e{e}.npz")
    return [np.load(f, allow_pickle=True) for f in sorted(glob.glob(pat))]


def emean(datas, key):
    return np.mean([d[key] for d in datas], axis=0)


def series(subdir, prefix):
    m2, spl, b22 = [], [], []
    for e in ECHECKS:
        ds = group(subdir, prefix, e)
        cnt, a_s, m2g = (emean(ds, "ax_counts"), emean(ds, "ax_a"),
                         emean(ds, "ax_m2_residue"))
        sel = (cnt > 30) & (a_s > 0)
        w = (a_s * cnt)[sel]
        m2.append(np.sum(m2g[sel] * w) / np.sum(w))
        phi, k2 = emean(ds, "phi_line"), ds[0]["kappa2"]
        ev = float(ds[0]["e"])
        use = (k2 > 0) & (k2 <= KMAX_FRAC * (int(ds[0]["n"]) // 3)
                          * np.exp(0.5 * ev))
        spl.append(np.sum(phi[0][use]) / np.sum(phi[2][use]) - 1.0)
        r = np.diag(emean(ds, "R"))
        b22.append(r[1] / r.sum() - 1 / 3)
    return np.array(m2), np.array(spl), np.array(b22)


m2_p, sp_p, b_p = series("n128", "")
m2_a, sp_a, b_a = series("n128axi", "axi_")

fig, (axl, axr) = plt.subplots(1, 2, figsize=(9.4, 3.7))

axl.axhline(0.133, color="0.7", lw=0.8, ls=":", zorder=0)
axl.text(0.02, 0.138, "isotropic floor", fontsize=7.5, color="0.5")
axl.plot(EVAL, m2_p, "o-", color=C_PLANE, lw=1.8, ms=5,
         label="plane strain  (A1 broken)")
axl.plot(EVAL, m2_a, "s-", color=C_AXI, lw=1.8, ms=5,
         label="axisymmetric  (A1 exact)")
axl.set_ylabel(r"azimuthal $m{=}2$ residue of $\Phi_{22}$", fontsize=9.5)
axl.set_title("Direct A1-violation measure", fontsize=10)
axl.legend(fontsize=8.5, frameon=False, loc="upper left")

axr.axhline(0.0, color="0.7", lw=0.8, ls=":", zorder=0)
axr.plot(EVAL, sp_p, "o-", color=C_PLANE, lw=1.8, ms=5,
         label=r"plane: $\phi_{11}/\phi_{33}{-}1$")
axr.plot(EVAL, sp_a, "s-", color=C_AXI, lw=1.8, ms=5,
         label=r"axisym: $\phi_{11}/\phi_{33}{-}1$")
axr.plot(EVAL, b_p, "^--", color=C_B22, lw=1.5, ms=5,
         label=r"$b_{22}$ (upwash) -- both geometries")
axr.plot(EVAL, b_a, "^:", color=C_B22, lw=1.5, ms=4, alpha=0.7)
axr.set_ylabel("transverse splitting  /  upwash anisotropy", fontsize=9.5)
axr.set_title("On-line falsifier + the robust quantity", fontsize=10)
axr.legend(fontsize=8, frameon=False, loc="lower left")

for ax in (axl, axr):
    ax.set_xlabel(r"accumulated strain $e = \int S\,dt$", fontsize=9.5)
    ax.tick_params(labelsize=8.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.93", lw=0.6, zorder=0)

fig.suptitle(r"The single-line (A1) closure error is a strain-geometry "
             r"effect, not an estimator error  ($Sk/\varepsilon=16$, "
             r"$128^3$, 4 seeds)", fontsize=10.5, y=1.02)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, "fig_a1_control." + ext),
                bbox_inches="tight", dpi=180)
print("saved fig_a1_control.pdf/.png")
