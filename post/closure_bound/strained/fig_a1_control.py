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
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                ".."))
from figstyle_jfm import FULL, panel, plt, save  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ECHECKS = ("0", "0.25", "0.5", "0.75", "1")
EVAL = [float(e) for e in ECHECKS]
KMAX_FRAC = 0.85
RATIO = "16"

# monochrome: geometry -> line style + marker (plane solid/circle,
# axisymmetric dashed/square); b22 (robust across geometry) grey dash-dot.
PLANE = dict(color="k", ls="-", marker="o", ms=3.4, mfc="k", mew=0.6)
AXI = dict(color="k", ls=(0, (5, 2)), marker="s", ms=3.4, mfc="none",
           mew=0.8)


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

fig, (axl, axr) = plt.subplots(1, 2, figsize=(FULL, 2.2))

axl.axhline(0.133, color="0.6", lw=0.6, ls=(0, (1, 1.2)), zorder=0)
axl.plot(EVAL, m2_p, lw=1.1, label="plane strain (A1 broken)", **PLANE)
axl.plot(EVAL, m2_a, lw=1.0, label="axisymmetric (A1 exact)", **AXI)
# label the floor in the empty band below the line (nothing lies below
# the floor), centred with clearance so it is not squeezed on the axis
axl.set_ylim(0.10, 0.41)
axl.text(0.32, 0.114, "isotropic floor", fontsize=7.5, color="0.4",
         va="center")
axl.set_ylabel(r"azimuthal $m{=}2$ residue of $\Phi_{22}$")
axl.legend(loc="upper left")

axr.axhline(0.0, color="0.6", lw=0.6, ls=(0, (1, 1.2)), zorder=0)
axr.plot(EVAL, sp_p, lw=1.1,
         label=r"plane: $\phi_{11}/\phi_{33}{-}1$", **PLANE)
axr.plot(EVAL, sp_a, lw=1.0,
         label=r"axisym: $\phi_{11}/\phi_{33}{-}1$", **AXI)
# b22 (upwash) is geometry-independent: both curves in one grey style so
# they merge, showing the acoustically relevant amplification is robust.
axr.plot(EVAL, b_p, color="0.45", ls=(0, (4, 1.4, 1, 1.4)), lw=0.9,
         marker="^", ms=3.0, mfc="none", mew=0.7,
         label=r"$b_{22}$ (upwash), both geometries")
axr.plot(EVAL, b_a, color="0.45", ls=(0, (4, 1.4, 1, 1.4)), lw=0.9,
         marker="^", ms=3.0, mfc="none", mew=0.7)
axr.set_ylabel("transverse splitting / upwash anisotropy")
axr.legend(loc="lower left")

for j, ax in enumerate((axl, axr)):
    ax.set_xlabel(r"accumulated strain $e = \int S\,dt$")
    panel(ax, "ab"[j], x=0.9)

fig.tight_layout(pad=0.4)
save(fig, os.path.join(HERE, "fig_a1_control"))
