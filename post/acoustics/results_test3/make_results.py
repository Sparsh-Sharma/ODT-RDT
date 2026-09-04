#!/usr/bin/env python3
r"""
Test-3 results: scale-resolved component anisotropy of the homogeneousStrain2
case (64 realizations, seeds 22..85 via shift 0..63), for the reply to Alan.

Produces fig_test3_anisotropy.{png,pdf} and results_table.txt.

Estimators: A_band = <E2>_band,rlz / <Eperp>_band,rlz (ratio of ensemble means;
the mean-of-ratios is biased high by heavy-tailed realization noise), with 95%
bootstrap CIs over realizations.  Reference is the ODT-internal isotropy fixed
point rho_iso = 1 (component equality; see scale_anisotropy.anisotropy_function).

Bands: low = energy-containing (k2 in [30,100]); high = 2.5-3 octaves below
(k2 in [300,800], inside the resolved range; the FFT interpolation floor sits
above k ~ 1200).  At e = 0 the IC has no modes above k ~ 402, so the high band
is floor-dominated there and excluded from interpretation.

Run: python3 make_results.py [dump_root]   (default ../../../data/homogeneousStrain2)
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import odt_io                                          # noqa: E402
import scale_anisotropy as sa                          # noqa: E402

STRAINS = [0.0, 1.0, 2.0, 3.0, 3.9]
BAND_LO = (30.0, 100.0)
BAND_HI = (300.0, 800.0)
NU = 2048
NBOOT = 2000


def band_means(root, e):
    dirs = sorted(glob.glob(os.path.join(root, "data", "data_*"))) or [root]
    rows = []
    for d in dirs:
        t = odt_io.load_fit_target(d, strain=e, Nu=NU)
        k2 = t["k2"]
        Ep = 0.5 * (t["E1"] + t["E3"])
        E2 = t["E2"]
        row = []
        for b in (BAND_LO, BAND_HI):
            m = (k2 >= b[0]) & (k2 <= b[1])
            row += [np.mean(E2[m]), np.mean(Ep[m])]
        R = [np.trapezoid(t[c], k2) for c in ("E1", "E2", "E3")]
        row.append(R[1] / sum(R))
        rows.append(row)
    return np.array(rows)


def ratio_ci(num, den, rng):
    n = num.size
    A = num.mean() / den.mean()
    idx = rng.integers(0, n, (NBOOT, n))
    Ab = num[idx].mean(axis=1) / den[idx].mean(axis=1)
    return A, np.percentile(Ab, 2.5), np.percentile(Ab, 97.5)


def main(root):
    rng = np.random.default_rng(0)
    table = []
    for e in STRAINS:
        r = band_means(root, e)
        lo = ratio_ci(r[:, 0], r[:, 1], rng)
        hi = ratio_ci(r[:, 2], r[:, 3], rng)
        table.append((e, r.shape[0], lo, hi, r[:, 4].mean()))

    lines = [
        f"homogeneousStrain2, {table[0][1]} realizations, Nu={NU}",
        f"A_band = <E2>/<Eperp> (ratio of ensemble means), 95% bootstrap CI",
        f"low band k2 in {BAND_LO}, high band k2 in {BAND_HI}",
        "",
        f"{'e':>5} {'nrlz':>5} {'A_low':>22} {'A_high':>22} {'u2^2/2kt':>9}",
    ]
    for e, n, lo, hi, f22 in table:
        lines.append(f"{e:5.1f} {n:5d}  {lo[0]:6.3f} [{lo[1]:5.3f},{lo[2]:5.3f}]"
                     f"   {hi[0]:6.3f} [{hi[1]:5.3f},{hi[2]:5.3f}] {f22:9.4f}")
    lines += [
        "",
        "level-0 (no-cascade) LRR target at e=4: u2^2/2kt = 0.6006; iso = 1/3.",
        "downscale transmission A_high/A_low: "
        + ", ".join(f"e={e:.1f}: {hi[0]/lo[0]:.2f}"
                    for e, n, lo, hi, f in table[1:]),
    ]
    txt = "\n".join(lines)
    print(txt)
    out = os.path.join(os.path.dirname(__file__), "results_table.txt")
    open(out, "w").write(txt + "\n")

    # figure ------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
    for e in STRAINS:
        r = sa.scale_anisotropy(root, strain=e, ensemble=True, Nu=NU,
                                nbin=24, kmax=1200.0, ref="equal")
        ax1.semilogx(r["k2"], r["A"], "o-", ms=3, label=f"$e$ = {e:.1f}")
    ax1.axhline(1.0, color="k", lw=0.8, ls=":")
    for b, c in ((BAND_LO, "0.85"), (BAND_HI, "0.92")):
        ax1.axvspan(*b, color=c, zorder=0)
    ax1.set_xlabel(r"$k_2$")
    ax1.set_ylabel(r"$A(k_2) = E_2/E_\perp$")
    ax1.set_title("scale-resolved component anisotropy (64 rlz)")
    ax1.legend(fontsize=8)

    es = [t[0] for t in table]
    for j, (lab, mk) in enumerate((("low band (30-100)", "o"),
                                   ("high band (300-800)", "s"))):
        A = [t[2 + j][0] for t in table]
        lo = [t[2 + j][0] - t[2 + j][1] for t in table]
        hi = [t[2 + j][2] - t[2 + j][0] for t in table]
        ax2.errorbar(es, A, yerr=[lo, hi], fmt=mk + "-", capsize=3, label=lab)
    ax2.axhline(1.0, color="k", lw=0.8, ls=":")
    ax2.set_xlabel(r"total strain $e$")
    ax2.set_ylabel(r"$A$ (band mean, 95% CI)")
    ax2.set_title("large- vs small-scale anisotropy")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    stem = os.path.join(os.path.dirname(__file__), "fig_test3_anisotropy")
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}", dpi=200)
    print(f"figure saved: {stem}.png/.pdf")


if __name__ == "__main__":
    default = os.path.join(os.path.dirname(__file__),
                           "..", "..", "..", "data", "homogeneousStrain2")
    main(sys.argv[1] if len(sys.argv) > 1 else default)
