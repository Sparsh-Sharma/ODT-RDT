"""JHTDB HIT null test for the LOS closure-bound estimator (SCOPING.md O5).

Consumes lines.npz from pull_lines.py (ensemble of velocity lines along y
from isotropic1024coarse) and answers, with error bars:

  T1  are the two transverse line spectra equal (phi_11 = phi_33)?
  T2  are the cross line spectra zero?
  T3  does the classical isotropic derivative relation hold?
  T4  do the LP bounds (raw and gamma=0) bracket the Crow value
      Pi = (2/5, -2/5, 0) k_t for the nominal plane strain, and what is
      the noise-inflated indeterminacy relative to the synthetic
      isotropic calibration?
  T5  false-positive metrics: b_ij anisotropy and the spectral
      transverse-splitting z-scores must be at the noise floor.

Conventions: line axis = x_2 (JHTDB y); components (u1,u2,u3)=(u_x,u_y,u_z);
integer wavenumbers (2*pi box); phi_mn(k2>0) are TWO-SIDED densities, so
R_mn = 2 * sum_{k2>0} phi_mn(k2), matching the estimator's half-axis
convention (PI_PREFACTOR includes the factor 2).
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir))
from axisym_estimator import Grid, _trapz_weights, bound_pi  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
KMAX_USE = 300          # resolved range for quantitative statements
N_K2_LP = 24            # log-binned k2 points fed to the LP
N_KPERP = 120


def line_cross_spectra(lines):
    """Per-line two-sided cross-spectral densities P_mn(k), k=0..N/2.

    lines: (L, N, 3) real. Returns P: (L, 6, N/2+1) for pairs
    [(1,1),(2,2),(3,3),(1,2),(1,3),(2,3)] (real parts).
    """
    lines = lines - lines.mean(axis=(0, 1), keepdims=True)  # global mean off
    f = np.fft.rfft(lines, axis=1) / lines.shape[1]
    pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
    p = np.empty((lines.shape[0], 6, f.shape[1]))
    for q, (m, n) in enumerate(pairs):
        cross = np.real(f[:, :, m] * np.conj(f[:, :, n]))
        # interior modes appear at +-k: two-sided density = one-sided/2 -> P
        p[:, q, :] = cross
    return p


def jackknife(vals):
    """Mean and jackknife SE over axis 0."""
    n = vals.shape[0]
    mean = vals.mean(axis=0)
    loo = (vals.sum(axis=0)[None, ...] - vals) / (n - 1)
    se = np.sqrt((n - 1) / n * ((loo - mean) ** 2).sum(axis=0))
    return mean, se


def log_bin(k, y_mean, y_se, n_bins, k_lo, k_hi):
    """Average (mean, se) into log-spaced k bins; se shrinks with bin count."""
    edges = np.geomspace(k_lo, k_hi, n_bins + 1)
    kc, ym, ys = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (k >= a) & (k < b)
        if sel.sum() == 0:
            continue
        kc.append(np.exp(np.mean(np.log(k[sel]))))
        ym.append(y_mean[..., sel].mean(axis=-1))
        ys.append(np.sqrt((y_se[..., sel] ** 2).mean(axis=-1) / sel.sum()))
    return np.array(kc), np.stack(ym, axis=-1), np.stack(ys, axis=-1)


def main():
    data = np.load(os.path.join(HERE, "lines.npz"))
    lines = data["lines"].astype(np.float64)
    n_lines = lines.shape[0]
    print(f"lines: {lines.shape}  dataset={data['dataset']}")

    p_all = line_cross_spectra(lines)
    k = np.arange(p_all.shape[2])          # integer wavenumbers
    mean, se = jackknife(p_all)            # (6, N/2+1)
    use = (k >= 1) & (k <= KMAX_USE)

    phi11, phi22, phi33 = mean[0], mean[1], mean[2]
    se11, se22, se33 = se[0], se[1], se[2]

    # --- component energies and anisotropy (T5) ---
    r_ii = 2.0 * mean[:3, 1:].sum(axis=1)
    kt = 0.5 * r_ii.sum()
    b_diag = r_ii / (2.0 * kt) - 1.0 / 3.0
    print(f"\nk_t = {kt:.4f}   R_ii = {np.round(r_ii, 4)}")
    print(f"b_11, b_22, b_33 = {np.round(b_diag, 4)}  (null target: ~0)")

    # --- T1: transverse splitting ---
    split = (phi11 - phi33)[use]
    split_se = np.sqrt(se11[use] ** 2 + se33[use] ** 2)
    z = split / split_se
    frac_2sig = float((np.abs(z) > 2.0).mean())
    print(f"\nT1 transverse splitting: mean|z| = {np.abs(z).mean():.2f}, "
          f"frac(|z|>2) = {frac_2sig:.3f}  (null target: ~0.05)")

    # --- T2: cross spectra (coherence) ---
    for q, name in ((3, "12"), (4, "13"), (5, "23")):
        coh = mean[q][use] / np.sqrt(phi11[use] * phi22[use])
        zq = mean[q][use] / se[q][use]
        print(f"T2 phi_{name}: mean coherence = {coh.mean():+.4f}, "
              f"mean|z| = {np.abs(zq).mean():.2f}")

    # --- T3: derivative relation (log-binned) ---
    kc, phim, _ = log_bin(k[use], np.stack([phi11[use], phi22[use],
                                            phi33[use]]),
                          np.stack([se11[use], se22[use], se33[use]]),
                          18, 1.0, KMAX_USE)
    dphi_l = np.gradient(phim[1], np.log(kc)) / kc
    pred_t = 0.5 * (phim[1] - kc * dphi_l)
    meas_t = 0.5 * (phim[0] + phim[2])
    rel = np.abs(meas_t - pred_t) / pred_t
    print(f"\nT3 derivative relation: median rel. deviation = "
          f"{np.median(rel):.3f} (interior bins "
          f"{np.median(rel[2:-2]):.3f})")

    # --- T4: estimator bounds ---
    kc_lp, phi_lp, se_lp = log_bin(
        k[use],
        np.stack([0.5 * (phi11[use] + phi33[use]), phi22[use]]),
        np.stack([0.5 * np.sqrt(se11[use] ** 2 + se33[use] ** 2),
                  se22[use]]),
        N_K2_LP, 1.0, KMAX_USE)
    band = float(np.median((se_lp / phi_lp).max(axis=0)) * 3.0)
    print(f"\nT4 LP: {kc_lp.size} k2 points, data band = {band:.3f}")

    kperp = np.geomspace(0.5, 2000.0, N_KPERP)
    grid = Grid(k2=kc_lp, kperp=kperp, w2=_trapz_weights(kc_lp),
                wperp=2.0 * np.pi * kperp * _trapz_weights(kperp))

    for label, kw in (("raw", {}),
                      ("gamma=0", dict(polarization_cap=0.0)),
                      ("gamma=0, lambda=4", dict(polarization_cap=0.0,
                                                 slope_cap=4.0))):
        bounds = bound_pi(grid, phi_lp[0], phi_lp[1], data_band=band, **kw)
        lo, hi = bounds.integrated(grid.w2)
        ok = bounds.status_ok.all()
        crow = np.array([0.4, -0.4, 0.0]) * kt
        bracket = bool(np.all(lo <= crow + 5e-3 * kt)
                       and np.all(crow - 5e-3 * kt <= hi))
        width = (hi - lo) / kt
        print(f"  [{label:18s}] solved={ok}  bracket_crow={bracket}  "
              f"width/kt = {np.round(width, 3)}")
        print(f"     lo/kt = {np.round(lo / kt, 3)}   "
              f"hi/kt = {np.round(hi / kt, 3)}")

    np.savez(os.path.join(HERE, "null_test_results.npz"),
             k=k, phi_mean=mean, phi_se=se, r_ii=r_ii, kt=kt,
             b_diag=b_diag, z_split=z, kc_lp=kc_lp, phi_lp=phi_lp,
             band=band, n_lines=n_lines)
    print("\nsaved null_test_results.npz")


if __name__ == "__main__":
    main()
