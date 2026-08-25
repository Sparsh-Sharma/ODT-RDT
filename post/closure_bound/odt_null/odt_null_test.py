"""ODT HIT null test: run the ODT nullHIT line ensemble (strain machinery on,
Astrain = 0, eddies on) through the same LOS estimator pipeline as the JHTDB
data — Alan's false-positive check on ODT itself.

Consumes nullhit_ensemble.npz (from extract_ensemble.py on caro).

Differences from the JHTDB pipeline, both forced by ODT's WALL BCs:
  * the line is not periodic and not homogeneous near the walls -> analyze
    the interior 80% segment with a Hann window (power-normalized);
  * wavenumbers are kappa_j = 2*pi*j / L_seg (physical units, L = 1 box).

Reports, per dump time: b_ij, T1 transverse-splitting z-scores, T2 cross
spectra; and at selected times the LP bounds (raw / gamma=0 / gamma=0+
lambda=4) with the Crow bracket check for the nominal plane strain.
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir))
sys.path.insert(0, os.path.join(HERE, os.pardir, "jhtdb_null"))
from axisym_estimator import Grid, _trapz_weights, bound_pi  # noqa: E402
from null_test import jackknife, log_bin  # noqa: E402

SEG_FRAC = 0.8          # interior segment (avoid wall-pinned ends)
N_K2_LP = 20
N_KPERP = 120
LP_TIMES = (2.0, 4.0)   # dump times at which to run the LP


def segment_spectra(lines, L):
    """Windowed cross-spectral densities on the interior segment.

    lines: (R, N, 3). Returns k (physical), P: (R, 6, M) two-sided densities
    such that R_mn ~= 2 * sum_{k>0} P * dk, dk = 2*pi/L_seg.
    """
    n = lines.shape[1]
    i0 = int(n * (1.0 - SEG_FRAC) / 2.0)
    seg = lines[:, i0:n - i0, :]
    nseg = seg.shape[1]
    l_seg = L * SEG_FRAC
    seg = seg - seg.mean(axis=(0, 1), keepdims=True)
    win = np.hanning(nseg)
    norm = np.sqrt((win ** 2).mean())
    f = np.fft.rfft(seg * win[None, :, None], axis=1) / (nseg * norm)
    pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
    dk = 2.0 * np.pi / l_seg
    p = np.empty((seg.shape[0], 6, f.shape[1]))
    for q, (m, mm) in enumerate(pairs):
        p[:, q, :] = np.real(f[:, :, m] * np.conj(f[:, :, mm])) / dk
    k = np.arange(f.shape[1]) * dk
    return k, p, dk


def main():
    data = np.load(os.path.join(HERE, "nullhit_ensemble.npz"))
    lines_all, times, L = data["lines"], data["times"], float(data["L"])
    print(f"ensemble: {lines_all.shape}, times {times}")

    kmax_use = 2.0 * np.pi / L * 60.0   # resolved: dxmax=0.02 -> j ~ 25; keep
    #  margin to show the tail but weight the LP toward resolved scales

    for di, t in enumerate(times):
        lines = lines_all[di]
        good = ~np.isnan(lines[:, 0, 0])
        lines = lines[good]
        k, p_all, dk = segment_spectra(lines.astype(np.float64), L)
        mean, se = jackknife(p_all)
        use = (k >= dk) & (k <= kmax_use)

        r_ii = 2.0 * mean[:3][:, k > 0].sum(axis=1) * dk
        kt = 0.5 * r_ii.sum()
        b = r_ii / (2.0 * kt) - 1.0 / 3.0
        z = ((mean[0] - mean[2]) / np.sqrt(se[0] ** 2 + se[2] ** 2))[use]
        zc = (mean[3:] / se[3:])[:, use]
        print(f"\n== t = {t:.1f}  (R={lines.shape[0]})  k_t = {kt:.4f}")
        print(f"   b_diag = {np.round(b, 4)}")
        print(f"   T1 split: mean|z| = {np.abs(z).mean():.2f}, "
              f"frac(|z|>2) = {(np.abs(z) > 2).mean():.3f}")
        print(f"   T2 cross: mean|z| 12/13/23 = "
              f"{np.round(np.abs(zc).mean(axis=1), 2)}")

        if not any(np.isclose(t, lt) for lt in LP_TIMES):
            continue

        kc_lp, phi_lp, se_lp = log_bin(
            k[use],
            np.stack([0.5 * (mean[0] + mean[2])[use], mean[1][use]]),
            np.stack([0.5 * np.sqrt(se[0] ** 2 + se[2] ** 2)[use],
                      se[1][use]]),
            N_K2_LP, dk, kmax_use)
        band = float(np.median((se_lp / phi_lp).max(axis=0)) * 3.0)
        kperp = np.geomspace(kc_lp[0] / 10.0, kc_lp[-1] * 10.0, N_KPERP)
        grid = Grid(k2=kc_lp, kperp=kperp, w2=_trapz_weights(kc_lp),
                    wperp=2.0 * np.pi * kperp * _trapz_weights(kperp))
        print(f"   LP ({kc_lp.size} pts, band {band:.3f}):")
        for label, kw in (("raw", {}),
                          ("gamma=0", dict(polarization_cap=0.0)),
                          ("g=0,l=4", dict(polarization_cap=0.0,
                                           slope_cap=4.0))):
            bd = bound_pi(grid, phi_lp[0], phi_lp[1], data_band=band, **kw)
            lo, hi = bd.integrated(grid.w2)
            crow = np.array([0.4, -0.4, 0.0]) * kt
            br = bool(np.all(lo <= crow + 5e-3 * kt)
                      and np.all(crow - 5e-3 * kt <= hi))
            print(f"     [{label:8s}] solved={bd.status_ok.all()} "
                  f"bracket_crow={br} lo/kt={np.round(lo / kt, 3)} "
                  f"hi/kt={np.round(hi / kt, 3)}")


if __name__ == "__main__":
    main()
