#!/usr/bin/env python3
r"""
Scale-resolved anisotropy diagnostic for homogeneous-strain ODT lines
(paper 2 / Test-3 discussion with Alan, 2026-09).

Purpose: quantify Alan's rate-competition diagnosis.  The triplet map compresses
scale by 3 per application (3^N after N maps) while transporting the component
imbalance down-scale unchanged; only the kernel events move energy between
components.  The observable consequence is how fast the *spectral* component
anisotropy decays with wavenumber -- this module measures exactly that, per
octave of k2 and per "cascade step" (one triplet map = a factor 3 in scale =
log2(3) octaves).

Method
------
The ODT line lies along x2, so E_2(k2) is the LONGITUDINAL 1-D spectrum and
E_1, E_3 are TRANSVERSE.  Even in perfect isotropy their ratio is not 1:

    rho_iso(k2) = E_2^iso / E_perp^iso   ->  2    (k2 -> 0)
                                          ->  3/4  (k2 -> inf, the 4/3 law),

so the raw ratio must be normalised by the isotropic reference before it means
"anisotropy".  We define

    rho(k2)  = E_2 / E_perp,          E_perp = (E_1 + E_3)/2,
    A(k2)    = rho(k2) / rho_iso(k2; L_iso),

with L_iso from a 2-parameter isotropic von Karman fit to E_perp (or given).
A(k2) = 1 at scales that are componentwise isotropic; |A - 1| is the scale-
resolved anisotropy amplitude.  The decay rate is a log-log fit over a band:

    |A - 1| ~ k2^s   =>  factor 2^s per octave,  factor 3^s per triplet map.

3^s is the number to hold against Alan's argument: it is the measured
multiplicative anisotropy reduction accompanying one 3x scale compression.

NOTE ON CONVENTIONS: everything here uses only the c0 = 0 (exact von Karman)
branch of axisym_family, so this diagnostic is untouched by the open h(mu) /
C-tensor decisions flagged there.
"""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

import axisym_family as af
from axisym_family import Params
import odt_io

# fit-grade inner grid (matches fit_family; ~1e-6 accuracy, ample here)
_HALF = 22.0
_N = 400

LOG2_3 = np.log2(3.0)  # octaves per triplet map (one map = 3x scale reduction)


# ----------------------------------------------------------------------
# Ratios
# ----------------------------------------------------------------------
def component_ratio(E1, E2, E3):
    """rho(k2) = E_2 / E_perp with E_perp = (E_1 + E_3)/2 (long./transverse)."""
    return np.asarray(E2, float) / (0.5 * (np.asarray(E1, float)
                                           + np.asarray(E3, float)))


def iso_reference_ratio(k2, L, half=_HALF, n=_N):
    """rho_iso(k2) for isotropic vK with wavenumber scale L (A0 cancels)."""
    p = Params(A0=1.0, c0=0.0, L2=L, Lperp=L)
    E1 = af.E_component(1, k2, p, half=half, n=n)   # transverse (= E3)
    E2 = af.E_component(2, k2, p, half=half, n=n)   # longitudinal
    return np.atleast_1d(E2) / np.atleast_1d(E1)


def fit_isotropic(k2, Eperp, half=_HALF, n=_N):
    """2-parameter (A0, L) isotropic vK fit to the transverse spectrum E_perp.

    Only L matters for the reference ratio (A0 cancels); it is fit in log
    space like fit_family.  Returns (A0, L).
    """
    k2 = np.asarray(k2, float)
    logE = np.log(np.clip(np.asarray(Eperp, float), 1e-300, None))

    def resid(theta):
        A0, L = np.exp(theta)
        p = Params(A0=A0, c0=0.0, L2=L, Lperp=L)
        Em = af.E_component(1, k2, p, half=half, n=n)
        return np.log(np.clip(Em, 1e-300, None)) - logE

    kpk = k2[np.argmax(k2 * np.asarray(Eperp))]     # crude energy wavenumber
    theta0 = np.log([max(np.max(Eperp), 1e-12), max(kpk, 1e-6)])
    sol = least_squares(resid, theta0, method="lm", xtol=1e-12, ftol=1e-12)
    A0, L = np.exp(sol.x)
    return float(A0), float(L)


# ----------------------------------------------------------------------
# Scale-resolved anisotropy and its decay rate
# ----------------------------------------------------------------------
def anisotropy_function(k2, E1, E2, E3, L_iso=None):
    """A(k2) = rho / rho_iso.  Returns a dict with A, rho, rho_iso, L_iso."""
    k2 = np.asarray(k2, float)
    rho = component_ratio(E1, E2, E3)
    A0_iso = None
    if L_iso is None:
        Eperp = 0.5 * (np.asarray(E1, float) + np.asarray(E3, float))
        A0_iso, L_iso = fit_isotropic(k2, Eperp)
    rho_iso = iso_reference_ratio(k2, L_iso)
    return {"k2": k2, "A": rho / rho_iso, "rho": rho, "rho_iso": rho_iso,
            "L_iso": L_iso, "A0_iso": A0_iso}


@dataclass
class DecayRate:
    slope: float                # d log2|A-1| / d log2 k2 over the band
    factor_per_octave: float    # 2**slope
    factor_per_map: float       # 3**slope  (per one triplet map = 3x in scale)
    npts: int                   # points used in the fit
    kband: tuple                # (kmin, kmax) actually used


def decay_per_octave(k2, A, kband=None, floor=1e-3):
    """Fit |A - 1| ~ k2^s over `kband` (default: all k2) where |A-1| > floor.

    Returns a DecayRate; slope < 0 means anisotropy decays toward small
    scales, slope ~ 0 means it persists (the Test-3 concern).  factor_per_map
    = 3**s is the multiplicative anisotropy reduction accompanying one
    triplet-map-sized (3x) scale compression.
    """
    k2 = np.asarray(k2, float)
    d = np.abs(np.asarray(A, float) - 1.0)
    m = d > floor
    if kband is not None:
        m &= (k2 >= kband[0]) & (k2 <= kband[1])
    n = int(np.count_nonzero(m))
    if n < 4:
        return DecayRate(np.nan, np.nan, np.nan, n,
                         kband if kband else (np.nan, np.nan))
    x = np.log2(k2[m])
    y = np.log2(d[m])
    s = float(np.polyfit(x, y, 1)[0])
    used = (float(k2[m].min()), float(k2[m].max()))
    return DecayRate(s, 2.0 ** s, 3.0 ** s, n, used)


# ----------------------------------------------------------------------
# Driver: dumps -> scale-resolved anisotropy report
# ----------------------------------------------------------------------
def scale_anisotropy(root, strain=None, smag=1.0, Nu=2048, nbin=48,
                     kmin=None, kmax=None, kband=None, L_iso=None,
                     floor=1e-3):
    """Full diagnostic on a run: read the dump at (nearest) `strain`, band-limit
    and log-bin the component spectra, form A(k2), and fit its decay rate.

    Returns a dict: k2, E1, E2, E3, rho, rho_iso, A, L_iso, decay (DecayRate),
    plus the dump metadata (e, t, L, file).
    """
    d = odt_io.load_fit_target(root, strain=strain, smag=smag, Nu=Nu)
    k2, E1, E2, E3 = d["k2"], d["E1"], d["E2"], d["E3"]
    sel = np.ones_like(k2, dtype=bool)
    if kmin is not None:
        sel &= k2 >= kmin
    if kmax is not None:
        sel &= k2 <= kmax
    k2, E1, E2, E3 = k2[sel], E1[sel], E2[sel], E3[sel]
    if nbin is not None:
        k2, (E1, E2, E3) = odt_io.log_bin(k2, [E1, E2, E3], nbin)

    an = anisotropy_function(k2, E1, E2, E3, L_iso=L_iso)
    decay = decay_per_octave(k2, an["A"], kband=kband, floor=floor)
    out = {"k2": k2, "E1": E1, "E2": E2, "E3": E3, "decay": decay,
           "e": d["e"], "t": d["t"], "L": d["L"], "file": d["file"]}
    out.update(an)
    return out


def report(res):
    """One-paragraph text summary of a scale_anisotropy() result."""
    dc = res["decay"]
    lines = [
        f"dump {res['file']}  (e = {res['e']:.4g}, domain L = {res['L']:.4g})",
        f"isotropic reference scale L_iso = {res['L_iso']:.5g}",
        f"|A-1| slope = {dc.slope:+.3f} per octave over "
        f"k2 in [{dc.kband[0]:.3g}, {dc.kband[1]:.3g}]  ({dc.npts} pts)",
        f"anisotropy factor per octave      : {dc.factor_per_octave:.3f}",
        f"anisotropy factor per triplet map : {dc.factor_per_map:.3f}"
        f"   (one map = 3x scale compression)",
    ]
    return "\n".join(lines)


# ----------------------------------------------------------------------
# Demo: controlled synthetic decay, recovered
# ----------------------------------------------------------------------
def _demo():
    L = af.ke_from_length(1.0)
    p_iso = Params(A0=1.0, c0=0.0, L2=L, Lperp=L)
    k2 = af.default_k2_grid(p_iso, half=16, n=120)[::3]
    E1 = af.E_component(1, k2, p_iso, half=_HALF, n=_N)
    E2i = af.E_component(2, k2, p_iso, half=_HALF, n=_N)
    d0, r, k0 = 0.4, 0.5, L                     # imposed: |A-1| ~ k^-0.5
    E2 = E2i * (1.0 + d0 * (k2 / k0) ** (-r))
    an = anisotropy_function(k2, E1, E2, E1)
    dc = decay_per_octave(k2, an["A"], kband=(2 * L, 50 * L))
    print("controlled synthetic: |A-1| ~ k^-0.5 imposed")
    print(f"recovered slope/octave = {dc.slope:+.4f}  (expect -0.5)")
    print(f"factor per octave      = {dc.factor_per_octave:.4f} "
          f"(expect {2 ** -r:.4f})")
    print(f"factor per triplet map = {dc.factor_per_map:.4f} "
          f"(expect {3 ** -r:.4f})")


if __name__ == "__main__":
    _demo()
