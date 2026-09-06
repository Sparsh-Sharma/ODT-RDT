#!/usr/bin/env python3
"""Refresh the ODT (ISO) rows of threeway.npz with the fresh same-binary
baseline ensemble fw_S2I (caro 2026-09-08, post-dilatation-fix binary,
identical precursor protocol), at e = 1.0 and 0.5.  The archived S2_ISO
ensemble predates the dilatation fix; band means agree within noise
except in the tail-sensitive finest band (+0.147 -> +0.076).  RDT and
DNS entries are untouched.  TYPESw entries are removed (different
binary; not referenced in the manuscript).  Writes threeway.npz in
place (backup threeway_archived.npz kept once).
Usage: python refresh_threeway.py <dir with fw_s2i_ensemble.npz>
"""
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, os.pardir, "odt_alloc"))
os.environ.setdefault("ALLOC_SLOW", "S2")
import threeway as T  # noqa: E402
from alloc_tests import centroid, line_spectra  # noqa: E402

QS = ("db22", "db11", "db33", "dsplit", "split", "b022", "R22", "R11",
      "R33")


def odt_ens(npz_path, di, f):
    d = np.load(npz_path)
    k0, (q11, _), (q22, _), (q33, _) = line_spectra(
        d["lines"][0], d["Ldump"][0], all_components=True)
    k1, (p11, _), (p22, _), (p33, _) = line_spectra(
        d["lines"][di], d["Ldump"][di], all_components=True)
    kref = centroid(k0, q11 + q22 + q33)
    return T.observables(k0, [q11, q22, q33], k1, [p11, p22, p33], f,
                         12.0 * kref)


def main():
    fwdir = sys.argv[1]
    src = os.path.join(HERE, "threeway.npz")
    bak = os.path.join(HERE, "threeway_archived.npz")
    if not os.path.exists(bak):
        shutil.copy(src, bak)
    old = dict(np.load(src, allow_pickle=True))
    out = {k: v for k, v in old.items() if "TYPESw" not in k}
    for e, di in ((1.0, 4), (0.5, 2), (0.0, 0)):
        f = np.exp(e / 2)
        obs = odt_ens(os.path.join(fwdir, "fw_s2i_ensemble.npz"), di, f)
        for q in QS:
            out[f"ISO_e{e}_{q}"] = obs[q]
    np.savez(src, **out)
    print("refreshed threeway.npz; fresh ISO db22 at e=1 (bands 2-8):")
    print("  " + " ".join(f"{x:+.3f}" for x in out["ISO_e1.0_db22"][1:]))


if __name__ == "__main__":
    main()
