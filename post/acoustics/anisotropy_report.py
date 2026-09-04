#!/usr/bin/env python3
r"""
Scale-resolved anisotropy report for a homogeneous-strain ODT run
(paper 2 / Test-3 rate-competition diagnostic; see scale_anisotropy.py).

Workflow matches spectrum_diagnostic2.py:
    Spyder: set DUMP_DIR below, press Run.
    CLI   : python3 anisotropy_report.py <dump_dir>

For each strain snapshot it forms the scale-resolved anisotropy
A(k2) = (E2/Eperp)/rho_iso(k2) and fits the decay |A-1| ~ k2^s, reporting the
anisotropy factor per octave (2^s) and per triplet map (3^s, one map = 3x
scale compression) -- the numbers for the discussion with Alan.  A slope
s ~ 0 with |A-1| well above zero is the "Test-3 signature": component
anisotropy transported down-scale faster than kernel events can relax it.

Outputs: a per-strain table on stdout and fig_scale_anisotropy.{png,pdf}
(saved next to the current working directory), plus an interactive window
when a display is available.
"""
import os
import sys

import numpy as np

import scale_anisotropy as sa
from scale_anisotropy import decay_per_octave

# ======================================================================
DUMP_DIR = r'C:\Users\shar_sp\Documents\ODT-post\data\homogeneousStrain2'
SMAG     = 1.0                    # strain magnitude S (e = S t)
STRAINS  = [0.0, 1.0, 2.0, 3.0, 4.0]
NUNIFORM = 2048                   # uniform resampling resolution for the FFT
NBIN     = 48                     # log-bins for the spectra (None = no binning)
KBAND    = None                   # (kmin, kmax) for the decay fit;
                                  # None -> auto (2*L_iso, 50*L_iso) per snapshot
FLOOR    = 1e-2                   # |A-1| below this is treated as isotropic
# ======================================================================


def analyze(root):
    rows = []
    for e in STRAINS:
        res = sa.scale_anisotropy(root, strain=e, smag=SMAG, Nu=NUNIFORM,
                                  nbin=NBIN, kband=KBAND, floor=FLOOR)
        if KBAND is None:                       # auto band from the fitted scale
            band = (2.0 * res["L_iso"], 50.0 * res["L_iso"])
            res["decay"] = decay_per_octave(res["k2"], res["A"],
                                            kband=band, floor=FLOOR)
        rows.append(res)
        print(sa.report(res))
        print("-" * 72)
    return rows


def summary_table(rows):
    print(f"{'e':>6} {'L_iso':>9} {'slope/oct':>10} {'per octave':>11} "
          f"{'per map':>9} {'npts':>5}")
    for r in rows:
        d = r["decay"]
        print(f"{r['e']:6.2f} {r['L_iso']:9.4g} {d.slope:10.3f} "
              f"{d.factor_per_octave:11.3f} {d.factor_per_map:9.3f} {d.npts:5d}")
    print("\nper map = multiplicative anisotropy change accompanying one "
          "triplet-map-sized (3x)\nscale compression; ~1.0 means anisotropy "
          "rides the cascade undiminished (Test 3).")


def make_figure(rows, stem="fig_scale_anisotropy"):
    try:
        import matplotlib
        if not os.environ.get("DISPLAY") and os.name != "nt":
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as ex:                     # HPC nodes without matplotlib
        print(f"(plotting skipped: {ex})")
        return None

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
    for r in rows:
        ax1.semilogx(r["k2"], r["A"], label=f"e = {r['e']:.1f}")
        d = np.abs(r["A"] - 1.0)
        m = d > FLOOR
        if np.any(m):
            ax2.loglog(r["k2"][m], d[m], ".", ms=4, label=f"e = {r['e']:.1f}")
    ax1.axhline(1.0, color="k", lw=0.8, ls=":")
    ax1.set_xlabel(r"$k_2$"); ax1.set_ylabel(r"$A(k_2)$")
    ax1.set_title("scale-resolved anisotropy (1 = isotropic)")
    ax1.legend(fontsize=8)
    ax2.set_xlabel(r"$k_2$"); ax2.set_ylabel(r"$|A-1|$")
    ax2.set_title("anisotropy amplitude and decay slope")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}", dpi=200)
    print(f"figure saved: {stem}.png / .pdf")
    try:
        plt.show()
    except Exception:
        pass
    return fig


def main(root=None):
    root = root or (sys.argv[1] if len(sys.argv) > 1 else DUMP_DIR)
    print(f"dump dir: {root}\n" + "=" * 72)
    rows = analyze(root)
    summary_table(rows)
    make_figure(rows)
    return rows


if __name__ == "__main__":
    main()
