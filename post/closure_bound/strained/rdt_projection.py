"""What exact linear RDT does to the PROJECTED (line) spectra, from the
Cauchy-RDT companions of the 128^3 pilot (same initial fields, evolved under
exact RDT, projected on the x2 line), against (i) the rigid-translation
hypothesis of the manuscript's Sec. 4.2.1 and (ii) the nonlinear DNS.

Rigid translation: phi_nn(k2, e) = J * phi_nn(k2 e^{-e/2}, 0) up to a
k2-independent amplitude. Shape ratio R_nn(k2) = [phi(k2,e)/int phi(.,e)] /
[phi(k2 e^{-e/2},0)/int phi(.,0)] * e^{-e/2}: equal to 1 at every k2 iff the
projected spectrum translates rigidly.

Splitting: phi11/phi33 - 1 versus k2 for exact RDT and for the DNS at the same
e -- decides whether the scale dependence of the splitting is linear.
"""
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "n128")
KREF = 3.0            # e=0 centroid of the total line spectrum (integer box units)
EDGES = np.geomspace(0.3, 12.0, 9)


def load_mean(pattern):
    files = sorted(glob.glob(os.path.join(D, pattern)))
    ph = np.mean([np.load(f, allow_pickle=True)["phi_line"] for f in files], axis=0)
    k = np.load(files[0], allow_pickle=True)["kappa2"]
    return k, ph, len(files)


def band(x, y, sel):
    out = []
    for a, b in zip(EDGES[:-1], EDGES[1:]):
        m = (x >= a) & (x < b) & sel
        out.append(np.mean(y[m]) if m.any() else np.nan)
    return np.array(out)


def band_ratio(x, num, den, sel):
    out = []
    for a, b in zip(EDGES[:-1], EDGES[1:]):
        m = (x >= a) & (x < b) & sel
        out.append(num[m].sum() / den[m].sum() - 1 if m.any() else np.nan)
    return np.array(out)


def main():
    k0, ph0, n0 = load_mean("chk_r0.8_s*_e0_rdt.npz")
    res = {}
    print(f"RDT companions: {n0} seeds; bands of k2(e)/k_c(0), k_c(0)={KREF}")
    print("bands: " + " ".join(f"[{a:.2f},{b:.2f})" for a, b in zip(EDGES[:-1], EDGES[1:])))
    for e in (0.5, 1.0):
        f = np.exp(e / 2)
        etag = "0.5" if e == 0.5 else "1"
        k1, ph1, _ = load_mean(f"chk_r0.8_s*_e{etag}_rdt.npz")
        kd, phd, nd = load_mean(f"chk_r0.8_s*_e{etag}.npz")
        sel = (k1 > 0) & (k1 <= 0.85 * 42 * f)
        x = k1 / KREF
        print(f"\n=== e={e}: shape ratio R_nn (1 = rigid translation) ===")
        for c, name in ((1, "phi22"), (0, "phi11"), (2, "phi33")):
            s1 = ph1[c][sel] / np.trapz(ph1[c][sel], k1[sel])
            p0i = np.interp(k1[sel] / f, k0, ph0[c])
            s0 = p0i / np.trapz(ph0[c][(k0 > 0) & (k0 <= 0.85 * 42)], k0[(k0 > 0) & (k0 <= 0.85 * 42)])
            R = s1 / (s0 / f)
            res[f"R_{name}_e{etag}"] = band(x[sel], R, np.ones(sel.sum(), bool))
            print(f"  {name:6s} " + " ".join(f"{v:6.3f}" for v in res[f"R_{name}_e{etag}"]))
        sr = band_ratio(x, ph1[0], ph1[2], sel)
        sd = band_ratio(kd / KREF, phd[0], phd[2], sel)
        res[f"split_rdt_e{etag}"] = sr
        res[f"split_dns_e{etag}"] = sd
        ir = ph1[0][sel].sum() / ph1[2][sel].sum() - 1
        idn = phd[0][sel].sum() / phd[2][sel].sum() - 1
        print(f"  splitting phi11/phi33-1, exact RDT : " + " ".join(f"{v:+6.3f}" for v in sr) + f"   integrated {ir:+.3f}")
        print(f"  splitting phi11/phi33-1, DNS 0.8   : " + " ".join(f"{v:+6.3f}" for v in sd) + f"   integrated {idn:+.3f}  ({nd} seeds)")
        # upwash spectrum: DNS vs exact RDT, band by band
        rr = band_ratio(x, np.interp(k1, kd, phd[1]), ph1[1], sel)
        res[f"dns_over_rdt_phi22_e{etag}"] = rr + 1
        print(f"  phi22 DNS / phi22 exact-RDT        : " + " ".join(f"{v+1:6.3f}" for v in rr))
    np.savez(os.path.join(HERE, "rdt_projection.npz"), edges=EDGES, **res)

    # ---- figure ----
    xc = np.sqrt(EDGES[:-1] * EDGES[1:])
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.2, 3.7))
    for name, col in (("phi22", "#2a78d6"), ("phi11", "#eb6834"), ("phi33", "#1baf7a")):
        a.plot(xc, res[f"R_{name}_e1"], "o-", color=col, ms=4, lw=1.5,
               label=r"$\phi_{%s}$, $e{=}1$" % name[3:])
        a.plot(xc, res[f"R_{name}_e0.5"], "o--", color=col, ms=3, lw=1.0, alpha=0.6,
               label=r"$\phi_{%s}$, $e{=}0.5$" % name[3:])
    a.axhline(1, color="k", lw=1.0, ls=":", label="rigid translation (ODT kinematics)")
    a.set_xscale("log")
    a.set_xlabel(r"$\kappa_2(e)/\kappa_c(0)$", fontsize=9.5)
    a.set_ylabel(r"shape ratio $R_{nn}$ (exact RDT / rigid)", fontsize=9.5)
    a.set_title("Exact linear RDT does not translate the line spectra rigidly", fontsize=10)
    a.legend(fontsize=7.2, frameon=False, ncol=2)
    b.plot(xc, res["split_dns_e1"], "k-", lw=2.2, label=r"DNS, $Sk/\varepsilon{=}0.8$, $e{=}1$")
    b.plot(xc, res["split_rdt_e1"], "-", color="#8e44ad", lw=2.0, label="exact RDT, $e{=}1$")
    b.plot(xc, res["split_dns_e0.5"], "k--", lw=1.2, label=r"DNS, $e{=}0.5$")
    b.plot(xc, res["split_rdt_e0.5"], "--", color="#8e44ad", lw=1.2, label="exact RDT, $e{=}0.5$")
    b.axhline(0, color="0.75", lw=0.8, ls=":")
    b.set_xscale("log")
    b.set_xlabel(r"$\kappa_2(e)/\kappa_c(0)$", fontsize=9.5)
    b.set_ylabel(r"$\phi_{11}/\phi_{33}-1$", fontsize=9.5)
    b.set_title("Scale dependence of the transverse splitting: linear or not?", fontsize=10)
    b.legend(fontsize=7.5, frameon=False)
    for ax in (a, b):
        ax.tick_params(labelsize=8.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="0.93", lw=0.6)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "fig_rdt_projection." + ext), bbox_inches="tight", dpi=180)
    print("saved fig_rdt_projection.pdf/.png")


if __name__ == "__main__":
    main()
