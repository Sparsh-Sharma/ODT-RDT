#!/usr/bin/env python3
r"""
Baseline vs Option A (anisotropy-gated eddy acceptance) comparison for the
Test-3 diagnostic.  Usage:

    python3 compare_optionA.py <baseline_root> <optionA_root>

Defaults to ../../../data/homogeneousStrain{2,2A}.  Works identically on the
local verification ensembles and the 1024-realization CARO campaign.

For each strain it reports the band anisotropies A = <E2>/<Eperp> (ratio of
ensemble means, 95% bootstrap CIs) in the energy-containing band and a band
2.5-3 octaves below, plus u2^2/2kt, for both cases side by side.  Option A
"working" = A_high pulled toward 1 (and A_high/A_low reduced) relative to
baseline at matched strain.
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


def case_stats(root):
    rng = np.random.default_rng(0)
    out = []
    for e in STRAINS:
        r = band_means(root, e)
        out.append((e, r.shape[0],
                    ratio_ci(r[:, 0], r[:, 1], rng),
                    ratio_ci(r[:, 2], r[:, 3], rng),
                    r[:, 4].mean()))
    return out


def fmt(ci):
    return f"{ci[0]:5.2f} [{ci[1]:4.2f},{ci[2]:4.2f}]"


def main(base_root, optA_root):
    sb = case_stats(base_root)
    sa_ = case_stats(optA_root)
    lines = [
        f"baseline: {base_root}  ({sb[0][1]} rlz)",
        f"option A: {optA_root}  ({sa_[0][1]} rlz)",
        f"bands: low k2 in {BAND_LO}, high k2 in {BAND_HI}; "
        f"A = <E2>/<Eperp>, 95% bootstrap CI",
        "",
        f"{'':>5} | {'baseline':^42} | {'option A':^42}",
        f"{'e':>5} | {'A_low':^17} {'A_high':^17} {'u2/2kt':>6} |"
        f" {'A_low':^17} {'A_high':^17} {'u2/2kt':>6}",
    ]
    for b, a in zip(sb, sa_):
        e = b[0]
        lines.append(f"{e:5.1f} | {fmt(b[2])} {fmt(b[3])} {b[4]:6.3f} |"
                     f" {fmt(a[2])} {fmt(a[3])} {a[4]:6.3f}")
    lines += ["", "downscale transmission A_high/A_low (e>0):"]
    for b, a in zip(sb[1:], sa_[1:]):
        lines.append(f"  e={b[0]:3.1f}:  baseline {b[3][0]/b[2][0]:5.2f}   "
                     f"option A {a[3][0]/a[2][0]:5.2f}")
    txt = "\n".join(lines)
    print(txt)
    here = os.path.dirname(__file__)
    open(os.path.join(here, "optionA_table.txt"), "w").write(txt + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
    for e, style in zip([2.0, 3.9], ["-", "--"]):
        for root, lab, col in ((base_root, "baseline", "C0"),
                               (optA_root, "option A", "C3")):
            r = sa.scale_anisotropy(root, strain=e, ensemble=True, Nu=NU,
                                    nbin=24, kmax=1200.0, ref="equal")
            ax1.semilogx(r["k2"], r["A"], style, color=col,
                         label=f"{lab}, e={e:.1f}")
    ax1.axhline(1.0, color="k", lw=0.8, ls=":")
    ax1.set_xlabel(r"$k_2$"); ax1.set_ylabel(r"$A(k_2)=E_2/E_\perp$")
    ax1.set_title("scale-resolved anisotropy")
    ax1.legend(fontsize=8)

    es = STRAINS
    for stats, lab, col in ((sb, "baseline", "C0"), (sa_, "option A", "C3")):
        for j, (bl, mk) in enumerate((("low", "o"), ("high", "s"))):
            A = [t[2 + j][0] for t in stats]
            lo = [t[2 + j][0] - t[2 + j][1] for t in stats]
            hi = [t[2 + j][2] - t[2 + j][0] for t in stats]
            ax2.errorbar(es, A, yerr=[lo, hi], fmt=mk + ("-" if j == 0 else "--"),
                         color=col, capsize=3, ms=4,
                         label=f"{lab}, {bl} band")
    ax2.axhline(1.0, color="k", lw=0.8, ls=":")
    ax2.set_xlabel(r"total strain $e$"); ax2.set_ylabel("A (band mean, 95% CI)")
    ax2.set_title("baseline vs option A")
    ax2.legend(fontsize=7)
    fig.tight_layout()
    stem = os.path.join(here, "fig_optionA_compare")
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}", dpi=200)
    print(f"figure saved: {stem}.png/.pdf")


if __name__ == "__main__":
    dd = os.path.join(os.path.dirname(__file__), "..", "..", "..", "data")
    base = sys.argv[1] if len(sys.argv) > 1 else os.path.join(dd, "homogeneousStrain2")
    opta = sys.argv[2] if len(sys.argv) > 2 else os.path.join(dd, "homogeneousStrain2A")
    main(base, opta)
