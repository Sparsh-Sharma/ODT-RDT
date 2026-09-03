"""Three-way comparison of the PROJECTED (line) spectra under plane strain at
Sk/eps = 0.8: ODT (S=2 precursor ensembles, 1024 realizations), exact linear
RDT (Cauchy companions of the 128^3 pilot) and the nonlinear DNS, on common
bands of kappa_2(e)/kappa_c(0) (kappa_c = e=0 centroid of each system's total
line spectrum; every system dilates by exp(e/2)).

Observable: the strain-induced change of the spectral component anisotropy
    b_nn(kappa_2) = phi_nn / (phi_11+phi_22+phi_33) - 1/3,
    Delta b_nn(kappa_2, e) = b_nn(kappa_2, e) - b_nn(kappa_2 e^{-e/2}, 0),
i.e. relative to the system's own e=0 state with every mode mapped to its
strained wavenumber. The e=0 reference matters: for real isotropic
turbulence the longitudinal (phi_22) and transverse (phi_11 = phi_33) line
spectra differ at every kappa_2 (b_22(kappa_2, 0) runs from +0.11 to -0.23
across the resolved range in the DNS), whereas the ODT line is
component-symmetric at e=0. Viscous decay, which is component-blind,
cancels in b_nn -- unlike in the shape ratio to a rigid translation (also
computed and printed, but confounded by decay for the two viscous systems).
Under a wavenumber-uniform rapid operator Delta b_nn(kappa_2) is flat;
exact RDT and the DNS are not. Note: sums over kappa_2 > 0 are NOT the
moment b_nn for the box data (the kappa_2 = 0 plane carries
component-dependent energy); the per-band values are what is compared.

Panels: b_22 (upwash), b_11 and b_33, and the transverse splitting
phi11/phi33 - 1, at e=1 (solid) and e=0.5 (dashed, RDT and DNS only).
This is the figure that replaces the strain-on/strain-off narrative of the
manuscript's Sec. 4.2 (Referee 3.1, 1.1, 2).
"""
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir, "odt_alloc"))
os.environ.setdefault("ALLOC_SLOW", "S2")
from alloc_tests import centroid, line_spectra, load  # noqa: E402

D = os.path.join(HERE, "n128")
EDGES = np.geomspace(0.3, 12.0, 9)
XC = np.sqrt(EDGES[:-1] * EDGES[1:])
KMAX_DNS = 0.85 * 42          # dealiased cutoff of the 128^3 box (integer units)
ODT_MODES = ("ISO", "TYPESw")  # ISO = the model; TYPESw shown faint


def bandmean(x, y):
    return np.array([np.mean(y[(x >= a) & (x < b)]) if ((x >= a) & (x < b)).any() else np.nan
                     for a, b in zip(EDGES[:-1], EDGES[1:])])


def bandfrac(x, num, den):
    """sum(num)/sum(den) per band (energy-weighted fraction)."""
    return np.array([num[(x >= a) & (x < b)].sum() / den[(x >= a) & (x < b)].sum()
                     if ((x >= a) & (x < b)).any() else np.nan
                     for a, b in zip(EDGES[:-1], EDGES[1:])])


def shape_ratio(k0, p0, k1, p1, f, kmax0):
    """R(k1) = [p1/int p1] / [p0(k1/f)/int p0] / f on the k1 grid, k1/f <= kmax0."""
    sel1 = (k1 > 0) & (k1 / f <= kmax0)
    sel0 = (k0 > 0) & (k0 <= kmax0)
    s1 = p1[sel1] / np.trapz(p1[sel1], k1[sel1])
    s0 = np.interp(k1[sel1] / f, k0, p0) / np.trapz(p0[sel0], k0[sel0])
    return k1[sel1], s1 / (s0 / f)


def observables(k0, P0, k1, P1, f, kmax0):
    """P = [phi11, phi22, phi33] arrays on k grids. Returns dict of band arrays."""
    kref = centroid(k0, P0[0] + P0[1] + P0[2])
    sel = (k1 > 0) & (k1 / f <= kmax0)
    x = k1[sel] / kref
    tot = P1[0][sel] + P1[1][sel] + P1[2][sel]
    # e=0 reference on the mapped grid: each mode at k1 came from k1/f
    P0i = [np.interp(k1[sel] / f, k0, P0[c]) for c in range(3)]
    tot0 = P0i[0] + P0i[1] + P0i[2]
    out = {"kref": kref}
    for c, name in ((1, "22"), (0, "11"), (2, "33")):
        b_e = bandfrac(x, P1[c][sel], tot) - 1.0 / 3.0
        b_0 = bandfrac(x, P0i[c], tot0) - 1.0 / 3.0
        out["b" + name] = b_e
        out["b0" + name] = b_0
        out["db" + name] = b_e - b_0
        kk, R = shape_ratio(k0, P0[c], k1, P1[c], f, kmax0)
        out["R" + name] = bandmean(kk / kref, R)
    out["split"] = bandfrac(x, P1[0][sel], P1[2][sel]) - 1.0
    out["split0"] = bandfrac(x, P0i[0], P0i[2]) - 1.0
    out["dsplit"] = out["split"] - out["split0"]
    return out


def dns_like(pattern_e0, pattern_e, f):
    f0 = sorted(glob.glob(os.path.join(D, pattern_e0)))
    f1 = sorted(glob.glob(os.path.join(D, pattern_e)))
    p0 = np.mean([np.load(p, allow_pickle=True)["phi_line"] for p in f0], axis=0)
    p1 = np.mean([np.load(p, allow_pickle=True)["phi_line"] for p in f1], axis=0)
    k0 = np.load(f0[0], allow_pickle=True)["kappa2"]
    k1 = np.load(f1[0], allow_pickle=True)["kappa2"]
    return observables(k0, p0[:3], k1, p1[:3], f, KMAX_DNS)


def odt_like(mode, di, f):
    d = load(f"S2_{mode}")
    k0, (q11, _), (q22, _), (q33, _) = line_spectra(d["lines"][0], d["Ldump"][0], all_components=True)
    k1, (p11, _), (p22, _), (p33, _) = line_spectra(d["lines"][di], d["Ldump"][di], all_components=True)
    kref = centroid(k0, q11 + q22 + q33)
    return observables(k0, [q11, q22, q33], k1, [p11, p22, p33], f, 12.0 * kref)


def main():
    res = {}
    for e, etag, di in ((1.0, "1", 4), (0.5, "0.5", 2), (0.0, "0", 0)):
        f = np.exp(e / 2)
        res[("RDT", e)] = dns_like("chk_r0.8_s*_e0_rdt.npz", f"chk_r0.8_s*_e{etag}_rdt.npz", f)
        res[("DNS", e)] = dns_like("chk_r0.8_s*_e0.npz", f"chk_r0.8_s*_e{etag}.npz", f)
        for m in ODT_MODES:
            res[(m, e)] = odt_like(m, di, f)
    print("bands of k2(e)/k_c(0): " + " ".join(f"[{a:.2f},{b:.2f})" for a, b in zip(EDGES[:-1], EDGES[1:])))
    QS = ("db22", "db11", "db33", "dsplit", "split", "b022", "R22", "R11", "R33")
    for e in (1.0, 0.5):
        print(f"\n===== e = {e} =====  (db = b(e) - b(e=0 mapped); b022 = e=0 reference of b22)")
        for q in QS:
            print(f"  {q}:")
            for s in ("RDT", "DNS") + ODT_MODES:
                print(f"    {s:7s} " + " ".join(f"{x:+6.3f}" for x in res[(s, e)][q]))
    np.savez(os.path.join(HERE, "threeway.npz"), edges=EDGES,
             **{f"{s}_e{e}_{q}": res[(s, e)][q] for (s, e) in res for q in QS})

    # ---- figure ----
    COL = {"RDT": "#8e44ad", "DNS": "k", "ISO": "#2a78d6", "TYPESw": "#eda100"}
    LBL = {"RDT": "exact linear RDT, projected", "DNS": r"DNS $128^3$, $Sk/\varepsilon{=}0.8$",
           "ISO": r"ODT, $Sk/\varepsilon{=}0.8$, isotropic kernel", "TYPESw": r"ODT, eddy types $p_3{=}0.5$"}
    LW = {"RDT": 2.0, "DNS": 2.3, "ISO": 1.8, "TYPESw": 1.0}
    AL = {"RDT": 1, "DNS": 1, "ISO": 1, "TYPESw": 0.5}
    MK = {"RDT": "D", "DNS": "o", "ISO": "s", "TYPESw": "^"}
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.9))
    panels = ((axs[0], "db22", r"strain-induced upwash anisotropy $\Delta b_{22}(\kappa_2)$", r"$\Delta b_{22}(\kappa_2)$"),
              (axs[1], "db11", r"$\Delta b_{11}(\kappa_2)$ (solid) and $\Delta b_{33}(\kappa_2)$ (dotted)", r"$\Delta b_{11},\ \Delta b_{33}$"),
              (axs[2], "split", r"transverse splitting (0 at $e{=}0$ for all)", r"$\phi_{11}/\phi_{33}-1$"))
    for ax, q, ttl, yl in panels:
        for s in ("RDT", "DNS", "ISO", "TYPESw"):
            ax.plot(XC, res[(s, 1.0)][q], "-", color=COL[s], lw=LW[s], alpha=AL[s], marker=MK[s], ms=3.8,
                    label=(LBL[s] + ", $e{=}1$") if ax is axs[0] else None)
            if s in ("RDT", "DNS"):
                ax.plot(XC, res[(s, 0.5)][q], "--", color=COL[s], lw=LW[s] * 0.55, alpha=0.8,
                        label=(LBL[s].split(",")[0] + ", $e{=}0.5$") if ax is axs[0] else None)
            if q == "db11":
                ax.plot(XC, res[(s, 1.0)]["db33"], ":", color=COL[s], lw=LW[s], alpha=AL[s], marker=MK[s], ms=3)
        ax.axhline(0, color="0.6", lw=0.9, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel(r"$\kappa_2(e)/\kappa_c(0)$", fontsize=9.5)
        ax.set_ylabel(yl, fontsize=9.5)
        ax.set_title(ttl, fontsize=9.8)
        ax.tick_params(labelsize=8.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="0.93", lw=0.6)
    axs[0].legend(fontsize=7, frameon=False, loc="lower left")
    fig.suptitle(r"Plane strain at $Sk/\varepsilon=0.8$: strain-induced anisotropy of the line spectra, "
                 "each relative to its own $e{=}0$ state -- ODT against exact projected RDT and DNS",
                 fontsize=10.5, y=1.02)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "fig_threeway." + ext), bbox_inches="tight", dpi=180)
    print("saved fig_threeway.pdf/.png")


if __name__ == "__main__":
    main()
