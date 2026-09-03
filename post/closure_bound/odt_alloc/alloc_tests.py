"""The three discriminating tests for the kernel-allocation modes
(notes/allocation_derivation.tex, Sec. 5).

Inputs: <case>_ensemble.npz from extract_ensemble.py for the 8 cases
  S8_{ISO,CHI,TYPESeq,TYPESw}   (rapid, S=8:   e=1 at t=0.125)
  S05_{ISO,CHI,TYPESeq,TYPESw}  (slow,  S=0.5: e=1 at t=2.0)
each (n_dumps=5 at e=0,.25,.5,.75,1; n_rlz; N; 3) with per-dump line
length Ldump (the line dilates under A_22).

Test 1  rapid limit: db22/de|0 and b22(e=1) per mode -> must be 2/15 and
        ~0.123 for EVERY mode (slow allocation must not touch the rapid
        response).
Test 2  moment level at S=0.5: b_ij(e) per mode vs the Level-0 table.
Test 3  spectral (decisive): phi_11(k2), phi_33(k2) and the splitting
        phi_11/phi_33 - 1 at S=0.5, e=1, per mode, against the 128^3 DNS
        line spectra at Sk/eps=0.8, e=1.  ISO predicted kappa-rigid;
        CHI/TYPES kappa-dependent.
"""

import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir, "jhtdb_null"))
sys.path.insert(0, os.path.join(HERE, os.pardir, "odt_null"))
from null_test import jackknife  # noqa: E402

MODES = ("ISO", "CHI", "TYPESeq", "TYPESw")
ECHK = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
# case-name prefixes of the rapid and slow ensembles; the 2026-09-03 campaign
# without precursor was S8/S05, the precursor campaign (tStrainOn=0.4) is
# S40/S2. OUT suffixes the result files.
RAPID = os.environ.get("ALLOC_RAPID", "S8")
SLOW = os.environ.get("ALLOC_SLOW", "S05")
OUT = os.environ.get("ALLOC_OUT", "")
SEG_FRAC = 0.8
DNS_DIR = os.path.join(HERE, os.pardir, "strained", "n128")


def load(case):
    f = os.path.join(HERE, f"{case.lower()}_ensemble.npz")
    return np.load(f) if os.path.exists(f) else None


def good_rlz(lines):
    """Realizations without a dump at this checkpoint are NaN-padded by the
    extractor (4-5 of 1024 at S=8); drop them."""
    return lines[~np.isnan(lines).any(axis=(1, 2))]


def b_of_dump(lines):
    """Component anisotropy with jackknife SE from a (R, N, 3) line set,
    fluctuations about the line mean, interior segment."""
    lines = good_rlz(lines)
    n = lines.shape[1]
    i0 = int(n * (1 - SEG_FRAC) / 2)
    seg = lines[:, i0:n - i0, :].astype(np.float64)
    seg = seg - seg.mean(axis=1, keepdims=True)
    r = (seg ** 2).mean(axis=1)                     # (R, 3) per-rlz R_ii
    nr = r.shape[0]

    def bf(rr):
        return rr / rr.sum(axis=-1)[..., None] - 1 / 3

    b = bf(r.mean(axis=0))
    loo = np.array([bf((r.sum(0) - r[i]) / (nr - 1)) for i in range(nr)])
    se = np.sqrt((nr - 1) / nr * ((loo - b) ** 2).sum(0))
    return b, se


def line_spectra(lines, ldump, all_components=False):
    """Windowed two-sided phi_11, phi_33 (per-rlz mean) vs physical k2.
    With all_components=True also returns phi_22."""
    lines = good_rlz(lines)
    n = lines.shape[1]
    i0 = int(n * (1 - SEG_FRAC) / 2)
    seg = lines[:, i0:n - i0, :].astype(np.float64)
    seg = seg - seg.mean(axis=1, keepdims=True)
    nseg = seg.shape[1]
    lseg = np.nanmean(ldump) * SEG_FRAC
    win = np.hanning(nseg)
    norm = np.sqrt((win ** 2).mean())
    f = np.fft.rfft(seg * win[None, :, None], axis=1) / (nseg * norm)
    dk = 2 * np.pi / lseg
    p11 = (np.abs(f[:, :, 0]) ** 2) / dk
    p33 = (np.abs(f[:, :, 2]) ** 2) / dk
    k = np.arange(f.shape[1]) * dk
    if all_components:
        p22 = (np.abs(f[:, :, 1]) ** 2) / dk
        return k, jackknife(p11), jackknife(p22), jackknife(p33)
    return k, jackknife(p11), jackknife(p33)


def centroid(k, p):
    """Energy-weighted spectral centroid over k>0 (the manuscript's own
    migration diagnostic); the common e=0 reference wavenumber for
    comparing ODT and DNS, whose absolute wavenumber units differ."""
    m = k > 0
    return float((k[m] * p[m]).sum() / p[m].sum())


def main():
    # ---------------- Test 1 + 2: moment level ----------------
    for tag, label in ((RAPID, f"Test 1 (rapid, {RAPID})"), (SLOW, f"Test 2 (slow, {SLOW})")):
        print(f"\n===== {label}: b_ij(e) per mode  [target rapid: slope 2/15=0.133, b22(1)~0.123]")
        print(f"{'mode':8s} " + "  ".join(f"e={e:<4}" + " " * 22 for e in ECHK))
        for mode in MODES:
            d = load(f"{tag}_{mode}")
            if d is None:
                print(f"{mode:8s} (missing)")
                continue
            row = []
            b22 = []
            for di in range(len(ECHK)):
                b, se = b_of_dump(d["lines"][di])
                b22.append(b[1])
                row.append(f"({b[0]:+.3f},{b[1]:+.3f},{b[2]:+.3f})")
            slope = (b22[1] - b22[0]) / 0.25
            print(f"{mode:8s} " + "  ".join(row) + f"   slope0={slope:.3f}")

    # ---------------- Test 3: spectral splitting at S=0.5, e=1 ----------------
    # Wavenumber convention: ODT and DNS are different nondimensional systems
    # (ODT IC built as a 1-D sum peaking at 8 waves per unit line; DNS IC a
    # 3-D shell spectrum with kp=4 per 2*pi box, whose 1-D line projection
    # peaks near 2). A 3-D peak is not available for ODT, so both are
    # normalized by the SAME measured 1-D quantity: the e=0 centroid of the
    # total line spectrum, x = kappa_2(e)/kappa_c(0). Both dilate by
    # exp(e/2), so the strain-induced shift is directly comparable.
    print(f"\n===== Test 3 (spectral, {SLOW}, e=1): phi_11/phi_33 - 1 in bands of "
          "k2/k_c(0)")
    bands = np.geomspace(0.3, 12.0, 9)
    results = {}
    dns1 = sorted(glob.glob(os.path.join(DNS_DIR, "chk_r0.8_s*_e1.npz")))
    dns0 = sorted(glob.glob(os.path.join(DNS_DIR, "chk_r0.8_s*_e0.npz")))
    if dns1 and dns0:
        ph = np.mean([np.load(f, allow_pickle=True)["phi_line"] for f in dns1], axis=0)
        ph0 = np.mean([np.load(f, allow_pickle=True)["phi_line"] for f in dns0], axis=0)
        k2d = np.load(dns1[0], allow_pickle=True)["kappa2"]
        k2d0 = np.load(dns0[0], allow_pickle=True)["kappa2"]
        kref_d = centroid(k2d0, ph0.sum(axis=0))
        use = (k2d > 0) & (k2d <= 0.85 * 42 * np.exp(0.5))
        xd = k2d / kref_d
        vals = []
        for a, b in zip(bands[:-1], bands[1:]):
            sel = (xd >= a) & (xd < b) & use
            vals.append(ph[0][sel].sum() / ph[2][sel].sum() - 1 if sel.any() else np.nan)
        results["DNS"] = np.array(vals)
        print(f"DNS k_c(0) = {kref_d:.3f} (integer box units)")
        print(f"{'DNS':8s} " + " ".join(f"{v:+.3f}" for v in vals)
              + f"   integrated={ph[0][use].sum()/ph[2][use].sum()-1:+.3f}")
    for mode in MODES:
        d = load(f"{SLOW}_{mode}")
        if d is None:
            continue
        k0, (q11, _), (q22, _), (q33, _) = line_spectra(
            d["lines"][0], d["Ldump"][0], all_components=True)
        kref = centroid(k0, q11 + q22 + q33)
        k, (p11, _), (p33, _) = line_spectra(d["lines"][-1], d["Ldump"][-1])
        x = k / kref
        vals = []
        for a, b in zip(bands[:-1], bands[1:]):
            sel = (x >= a) & (x < b)
            vals.append(p11[sel].sum() / p33[sel].sum() - 1 if sel.any() else np.nan)
        results[mode] = np.array(vals)
        print(f"{mode:8s} " + " ".join(f"{v:+.3f}" for v in vals)
              + f"   integrated={p11[x>0].sum()/p33[x>0].sum()-1:+.3f}"
              + f"   k_c(0)={kref/(2*np.pi):.2f} waves/L")
    print("bands (k2/k_c(0)): " + " ".join(f"[{a:.2f},{b:.2f})" for a, b in zip(bands[:-1], bands[1:])))
    np.savez(os.path.join(HERE, f"alloc_tests_results{OUT}.npz"), bands=bands,
             **{f"split_{m}": v for m, v in results.items()})


if __name__ == "__main__":
    main()
