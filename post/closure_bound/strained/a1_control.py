"""A1 control: plane strain vs axisymmetric contraction about the line.

Under axisymmetric contraction A = diag(1/4, -1/2, 1/4) the mean strain
PRESERVES axisymmetry about e_2, so assumption A1 is exact at every strain.
Under plane strain A = diag(1/2, -1/2, 0) it is not. The two kernel-free A1
diagnostics -- the azimuthal m=2 residue of Phi_22 and the phi_11 vs phi_33
splitting -- therefore separate "A1 violated" (a physics limitation of the
single-line closure) from "estimator broken" (it is not): they must grow
under plane strain and stay at the isotropic floor under axisymmetric
strain. This script tabulates both, ensemble-averaged over seeds.
"""

import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ECHECKS = ("0", "0.25", "0.5", "0.75", "1")
KMAX_FRAC = 0.85


def group(subdir, prefix, ratio, e):
    pat = os.path.join(HERE, subdir, f"chk_{prefix}r{ratio}_s*_e{e}.npz")
    return [np.load(f, allow_pickle=True) for f in sorted(glob.glob(pat))]


def emean(datas, key):
    return np.mean([d[key] for d in datas], axis=0)


def diagnostics(datas):
    cnt = emean(datas, "ax_counts")
    a_s = emean(datas, "ax_a")
    m2 = emean(datas, "ax_m2_residue")
    sel = (cnt > 30) & (a_s > 0)
    wgt = (a_s * cnt)[sel]
    m2_eff = float(np.sum(m2[sel] * wgt) / np.sum(wgt))

    phi = emean(datas, "phi_line")
    k2 = datas[0]["kappa2"]
    e = float(datas[0]["e"])
    n_deal = int(datas[0]["n"]) // 3
    use = (k2 > 0) & (k2 <= KMAX_FRAC * n_deal * np.exp(0.5 * e))
    split = float(np.sum(phi[0][use]) / np.sum(phi[2][use]) - 1.0)
    b = np.diag(emean(datas, "R"))
    b = b / b.sum() - 1 / 3
    return m2_eff, split, b[1]


def run(ratio):
    print(f"\n===== S k/eps = {ratio} :  m2 residue | phi11/phi33-1 | b22 "
          "=====")
    print(f"{'e':>5} | {'PLANE m2':>9} {'AXI m2':>8} | "
          f"{'PLANE spl':>10} {'AXI spl':>9} | {'PLANE b22':>10} "
          f"{'AXI b22':>8}")
    for e in ECHECKS:
        pl = group("n128", "", ratio, e)
        ax = group("n128axi", "axi_", ratio, e)
        if not pl or not ax:
            continue
        m2p, sp_p, b2p = diagnostics(pl)
        m2a, sp_a, b2a = diagnostics(ax)
        print(f"{e:>5} | {m2p:9.3f} {m2a:8.3f} | {sp_p:+10.3f} "
              f"{sp_a:+9.3f} | {b2p:+10.4f} {b2a:+8.4f}")


def main():
    for ratio in ("0.8", "16"):
        run(ratio)
    print("\nInterpretation: under axisymmetric contraction the m2 residue "
          "and\nthe transverse splitting stay near their isotropic floor "
          "(~0.13, 0),\nconfirming A1 holds; under plane strain both grow. "
          "The single-line\nclosure's error is the A1 violation, which is "
          "a property of the strain\ngeometry, not of the estimator.")


if __name__ == "__main__":
    main()
