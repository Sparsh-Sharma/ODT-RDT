"""First strained-turbulence pass through the closure-bound estimator.

Consumes the 128^3 pilot checkpoints (n128/chk_r{R}_s{S}_e{E}[_rdt].npz).
Per (ratio, e), ensemble-averaged over seeds:

  1. LP bounds [Pi-, Pi+] from the DNS line spectra alone (raw / gamma=0 /
     gamma=0+lambda=4) vs the EXACT Pi^(r) from the full spectral field --
     does the bound bracket the truth under real (non-axisymmetric) strain,
     and how wide is it?
  2. Measured polarization gamma_eff = |c|/a and the m=2 azimuthal residue
     of Phi_22 vs e -- adjudicates the stage-b polarization cap and the
     O(b^2) A1-error claim directly.
  3. On-line A1 diagnostics under strain: phi_11 vs phi_33 splitting.
  4. DNS vs the closed-form RDT companion: b_22(e).
"""

import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir))
from axisym_estimator import Grid, _trapz_weights, bound_pi  # noqa: E402

RATIOS = ("0.8", "16")
ECHECKS = ("0", "0.25", "0.5", "0.75", "1")
N_KPERP = 120
KMAX_FRAC = 0.85          # use k2 up to this fraction of the dealiased max


def load_group(ratio, e, rdt=False):
    suff = "_rdt" if rdt else ""
    pat = os.path.join(HERE, "n128", f"chk_r{ratio}_s*_e{e}{suff}.npz")
    files = sorted(glob.glob(pat))
    return [np.load(f, allow_pickle=True) for f in files]


def ens_mean(datas, key):
    return np.mean([d[key] for d in datas], axis=0)


def analyze(ratio):
    print(f"\n================ S k/eps = {ratio} ================")
    print("e     b22_DNS  b22_RDT   gamma_eff  m2_res   "
          "phi11/phi33-1")
    rows = {}
    for e in ECHECKS:
        dns = load_group(ratio, e)
        rdt = load_group(ratio, e, rdt=True)
        if not dns:
            continue
        r_mean = ens_mean(dns, "R")
        b_dns = np.diag(r_mean) / np.trace(r_mean) - 1 / 3
        r_rdt = ens_mean(rdt, "R")
        b_rdt = np.diag(r_rdt) / np.trace(r_rdt) - 1 / 3

        # polarization content and A1 residue from the (a,c) grids,
        # energy-weighted over well-populated bins
        cnt = ens_mean(dns, "ax_counts")
        a_s = ens_mean(dns, "ax_a")
        c_s = ens_mean(dns, "ax_c")
        m2 = ens_mean(dns, "ax_m2_residue")
        sel = (cnt > 30) & (a_s > 0)
        wgt = (a_s * cnt)[sel]
        # polarization ratio gamma_eff = <|c|>/<a>, energy weighted
        gamma_eff = float(np.sum(np.abs(c_s[sel]) * cnt[sel])
                          / np.sum(a_s[sel] * cnt[sel]))
        m2_eff = float(np.sum(m2[sel] * wgt) / np.sum(wgt))

        # on-line transverse splitting (A1 falsifier T1 under strain)
        phi = ens_mean(dns, "phi_line")          # (6, nk2)
        k2 = dns[0]["kappa2"]
        n_dealias = int(dns[0]["n"]) // 3
        use = (k2 > 0) & (k2 <= KMAX_FRAC * n_dealias
                          * np.exp(0.5 * float(e)))
        split = float(np.sum(phi[0][use]) / np.sum(phi[2][use]) - 1.0)

        rows[e] = dict(b_dns=b_dns, b_rdt=b_rdt, gamma=gamma_eff,
                       m2=m2_eff, split=split, phi=phi, k2=k2, use=use,
                       kt=float(ens_mean(dns, "kt")),
                       smag=float(ens_mean(dns, "smag")),
                       pi=ens_mean(dns, "pi_rapid"))
        print(f"{e:>4}  {b_dns[1]:+.4f}  {b_rdt[1]:+.4f}   "
              f"{gamma_eff:8.3f}  {m2_eff:.3f}    {split:+.3f}")

    # LP bounds vs exact Pi at selected strains
    print("\n  LP bounds vs exact Pi (per k_t S, S=1 units: values/kt):")
    for e in ("0.25", "0.5", "1"):
        if e not in rows:
            continue
        r = rows[e]
        kt = r["kt"]
        phi11 = 0.5 * (r["phi"][0] + r["phi"][2])[r["use"]]
        phi22 = r["phi"][1][r["use"]]
        k2v = r["k2"][r["use"]]
        # thin to ~24 log-spaced points
        idx = np.unique(np.geomspace(1, k2v.size - 1, 24).astype(int))
        grid = Grid(k2=k2v[idx], kperp=np.geomspace(0.3, 3 * k2v.max(),
                                                    N_KPERP),
                    w2=_trapz_weights(k2v[idx]),
                    wperp=None)
        grid.wperp = (2 * np.pi * grid.kperp
                      * _trapz_weights(grid.kperp))
        # estimator kernels are unit-S (A = diag(1/2,-1/2,0)); the DNS
        # Pi carries the actual S -> compare per unit k_t AND per unit S
        pi_ex = np.diag(r["pi"]) / (kt * r["smag"])
        print(f"  e={e} (S={r['smag']:.3f}):  exact Pi/(kt S) = "
              f"{np.round(pi_ex, 3)}")
        for label, kw in (("raw", {}),
                          ("gamma=0", dict(polarization_cap=0.0)),
                          ("g0,l4", dict(polarization_cap=0.0,
                                         slope_cap=4.0))):
            bd = bound_pi(grid, phi11[idx], phi22[idx], data_band=0.05,
                          **kw)
            lo, hi = bd.integrated(grid.w2)
            lo, hi = lo / kt, hi / kt
            inside = np.all((lo - 5e-3 <= pi_ex) & (pi_ex <= hi + 5e-3))
            print(f"    [{label:8s}] ok={bd.status_ok.all()} "
                  f"bracket={inside}  lo={np.round(lo, 3)} "
                  f"hi={np.round(hi, 3)}")
    return rows


def main():
    for ratio in RATIOS:
        analyze(ratio)


if __name__ == "__main__":
    main()
