"""b_ij(t) with jackknife errors for the three whitening variants, plus the
full null suite on the symmetric-whitening ensemble. Saves whitening_bt.npz
for the note's figure."""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir))
sys.path.insert(0, os.path.join(HERE, os.pardir, "jhtdb_null"))
from odt_null_test import segment_spectra  # noqa: E402

CASES = [("forward", "nullhit_ensemble.npz"),
         ("reversed", "nullhitrev_ensemble.npz"),
         ("symmetric", "nullhitsym_ensemble.npz")]


def b_of_t(fname):
    d = np.load(os.path.join(HERE, fname))
    lines_all, times, L = d["lines"], d["times"], float(d["L"])
    bs, ses = [], []
    for di in range(len(times)):
        lines = lines_all[di].astype(np.float64)
        k, p, dk = segment_spectra(lines, L)
        r = 2.0 * p[:, :3, 1:].sum(axis=2) * dk
        n = r.shape[0]

        def bfun(rr):
            return rr / rr.sum(axis=-1)[..., None] - 1.0 / 3.0

        b_full = bfun(r.mean(axis=0))
        loo = np.array([bfun((r.sum(axis=0) - r[i]) / (n - 1))
                        for i in range(n)])
        se = np.sqrt((n - 1) / n * ((loo - b_full) ** 2).sum(axis=0))
        bs.append(b_full)
        ses.append(se)
    return times, np.array(bs), np.array(ses)


def main():
    out = {}
    for label, fname in CASES:
        times, b, se = b_of_t(fname)
        out[label + "_b"], out[label + "_se"] = b, se
        out["times"] = times
        print(f"--- {label}")
        for di in (1, 4, 8):
            print(f"  t={times[di]:.1f}  " + "  ".join(
                f"b{c+1}{c+1}={b[di, c]:+.4f}+-{se[di, c]:.4f}"
                for c in range(3)))
    np.savez(os.path.join(HERE, "whitening_bt.npz"), **out)

    # symmetric ensemble: full spectral nulls at t = 0.5 and 2
    d = np.load(os.path.join(HERE, "nullhitsym_ensemble.npz"))
    lines_all, times = d["lines"], d["times"]
    from null_test import jackknife  # noqa: E402
    for di in (1, 4):
        lines = lines_all[di].astype(np.float64)
        k, p, dk = segment_spectra(lines, float(d["L"]))
        mean, se = jackknife(p)
        use = (k >= dk) & (k <= 2 * np.pi * 60)
        z = ((mean[0] - mean[2]) / np.sqrt(se[0] ** 2 + se[2] ** 2))[use]
        zc = (mean[3:] / se[3:])[:, use]
        print(f"sym t={times[di]:.1f}: T1 frac(|z|>2)={np.mean(np.abs(z) > 2):.3f}"
              f"  T2 mean|z| = {np.round(np.abs(zc).mean(axis=1), 2)}")
    print("saved whitening_bt.npz")


if __name__ == "__main__":
    main()
