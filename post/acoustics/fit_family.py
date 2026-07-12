#!/usr/bin/env python3
r"""
Fit the axisymmetric vK-stretch spectrum family (axisym_family.Params) to
component line spectra E_i(k2), for paper 2 (Gate A, Session 2).

The forward observables are the 1-D line spectra along the ODT-resolved axis,

    E_1(k2) = E_3(k2)   (perpendicular components; equal by axisymmetry)
    E_2(k2)             (the ODT/upwash-normal component)

which is exactly what a homogeneous-strain ODT run delivers after the standard
line-FFT reduction (cf. spectrum_diagnostic2.py: E_i(k) per velocity
component).  The four parameters (A0, c0, L2, Lperp) are recovered by a
log-space least-squares fit of the model E_1, E_2 to the target curves.

This module is DATA-INDEPENDENT: it does not read run data.  Point
`fit_family` at target arrays produced either synthetically (see
test_fit_family.py) or by the ODT reduction once those dumps are reachable.

    # DECISION NEEDED (inherited from axisym_family): the c0 != 0 branch uses
    # the P2-Legendre h(mu) stub and the rank-1 projected-axis C-tensor.  Fits
    # that lean on c0 (component-energy anisotropy) will shift if note sec.3/
    # paper1 sec.5 pin different forms.  L2 = Lperp isotropy and the A0/length
    # fit are unaffected by that choice.

Conventions and the spectral model live in axisym_family; see its docstring.
"""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

import axisym_family as af
from axisym_family import Params

# Cheaper log-grid for the inner E-integral during fitting (accuracy ~1e-6,
# ample for a least-squares target); the science-grade default lives in
# axisym_family (N_LOG=800).
_FIT_HALF = 22.0
_FIT_N = 400


# ----------------------------------------------------------------------
# Forward observables
# ----------------------------------------------------------------------
def line_spectra(p, k2_grid, half=_FIT_HALF, n=_FIT_N):
    """Model line spectra (E1, E2) on k2_grid.  E3 == E1 by axisymmetry."""
    E1 = af.E_component(1, k2_grid, p, half=half, n=n)
    E2 = af.E_component(2, k2_grid, p, half=half, n=n)
    return np.atleast_1d(E1), np.atleast_1d(E2)


def observables(p, k2_grid):
    """Full forward map params -> observables used downstream.

    Returns a dict with the line spectra (E1=E3, E2), the component energies
    R = [R11, R22, R33], and the turbulent kinetic energy k_t.
    """
    E1, E2 = line_spectra(p, k2_grid)
    R = af.R_ii(k2_grid, p)
    return {
        "E1": E1, "E2": E2, "E3": E1,
        "R": R,
        "kt": 0.5 * np.sum(R),
    }


# ----------------------------------------------------------------------
# Parameter <-> optimizer-vector packing (positivity via log; c0 free)
# ----------------------------------------------------------------------
def _pack(p):
    return np.array([np.log(p.A0), p.c0, np.log(p.L2), np.log(p.Lperp)])


def _unpack(theta):
    lnA0, c0, lnL2, lnLperp = theta
    return Params(A0=np.exp(lnA0), c0=c0, L2=np.exp(lnL2), Lperp=np.exp(lnLperp))


# ----------------------------------------------------------------------
# Residuals and fit
# ----------------------------------------------------------------------
def _residuals(theta, k2_grid, logE1_t, logE2_t, mask1, mask2, half, n):
    p = _unpack(theta)
    E1, E2 = line_spectra(p, k2_grid, half=half, n=n)
    # log-space residuals weight all decades of the vK spectrum equally.
    r1 = np.log(np.clip(E1[mask1], 1e-300, None)) - logE1_t[mask1]
    r2 = np.log(np.clip(E2[mask2], 1e-300, None)) - logE2_t[mask2]
    return np.concatenate([r1, r2])


def initial_guess(k2_grid, E1_t, E2_t):
    """Data-driven starting Params: energy scale from the perpendicular peak,
    amplitude from the integral, isotropic (c0=0, L2=Lperp) to start.
    """
    k2 = np.asarray(k2_grid, float)
    E1 = np.asarray(E1_t, float)
    E2 = np.asarray(E2_t, float)
    # peak of E2 (longitudinal) sits near the energy wavenumber
    kpk = k2[np.argmax(E2)] if np.any(E2 > 0) else k2[len(k2) // 2]
    L0 = max(kpk, 1e-6)
    # crude amplitude: match E2 peak height of the isotropic model at L0
    p_try = Params(A0=1.0, c0=0.0, L2=L0, Lperp=L0)
    _, E2m = line_spectra(p_try, k2)
    scale = np.max(E2) / max(np.max(E2m), 1e-300)
    return Params(A0=scale, c0=0.0, L2=L0, Lperp=L0)


@dataclass
class FitResult:
    params: Params
    success: bool
    cost: float
    nfev: int
    message: str


def fit_family(k2_grid, E1_target, E2_target, p0=None,
               rel_floor=1e-4, half=_FIT_HALF, n=_FIT_N,
               c0_bounds=(-2.0, 2.0)):
    """Fit (A0, c0, L2, Lperp) to target line spectra E1(=E3), E2.

    Parameters
    ----------
    k2_grid : positive log-spaced k2 nodes (see axisym_family.default_k2_grid).
    E1_target, E2_target : target line spectra on k2_grid (E1 = E3).
    p0 : optional starting Params (else data-driven initial_guess).
    rel_floor : points with E_target < rel_floor * max(E_target) are dropped
        (interpolation/round-off floor, as in spectrum_diagnostic2.py).

    Returns
    -------
    FitResult with the recovered Params and optimizer diagnostics.
    """
    k2_grid = np.asarray(k2_grid, float)
    E1_target = np.asarray(E1_target, float)
    E2_target = np.asarray(E2_target, float)

    mask1 = E1_target > rel_floor * E1_target.max()
    mask2 = E2_target > rel_floor * E2_target.max()
    logE1_t = np.log(np.clip(E1_target, 1e-300, None))
    logE2_t = np.log(np.clip(E2_target, 1e-300, None))

    if p0 is None:
        p0 = initial_guess(k2_grid, E1_target, E2_target)
    theta0 = _pack(p0)

    lo = np.array([-np.inf, c0_bounds[0], -np.inf, -np.inf])
    hi = np.array([np.inf, c0_bounds[1], np.inf, np.inf])

    sol = least_squares(
        _residuals, theta0, bounds=(lo, hi),
        args=(k2_grid, logE1_t, logE2_t, mask1, mask2, half, n),
        method="trf", xtol=1e-12, ftol=1e-12, max_nfev=400,
    )
    return FitResult(
        params=_unpack(sol.x),
        success=bool(sol.success),
        cost=float(sol.cost),
        nfev=int(sol.nfev),
        message=str(sol.message),
    )


# ----------------------------------------------------------------------
# Convenience: synthetic target generator (used by tests / demos)
# ----------------------------------------------------------------------
def synthetic_targets(p_true, k2_grid, noise=0.0, seed=0):
    """Generate (E1, E2) targets from known params, with optional lognormal
    multiplicative noise of relative level `noise`."""
    E1, E2 = line_spectra(p_true, k2_grid)
    if noise > 0.0:
        rng = np.random.default_rng(seed)
        E1 = E1 * np.exp(rng.normal(0.0, noise, size=E1.shape))
        E2 = E2 * np.exp(rng.normal(0.0, noise, size=E2.shape))
    return E1, E2


def _demo():
    p_true = Params(A0=1.4, c0=0.6, L2=0.8, Lperp=1.3)
    k2 = af.default_k2_grid(p_true, half=18, n=140)
    # thin to a realistic number of spectral points
    k2 = k2[::4]
    E1, E2 = synthetic_targets(p_true, k2, noise=0.0)
    res = fit_family(k2, E1, E2)
    pt, pf = p_true, res.params
    print("synthetic fit (noise-free):", res.message, "| nfev", res.nfev)
    print(f"{'param':>7} {'true':>10} {'fit':>10} {'rel.err':>10}")
    for name, a, b in [("A0", pt.A0, pf.A0), ("c0", pt.c0, pf.c0),
                       ("L2", pt.L2, pf.L2), ("Lperp", pt.Lperp, pf.Lperp)]:
        denom = abs(a) if abs(a) > 1e-9 else 1.0
        print(f"{name:>7} {a:10.5f} {b:10.5f} {abs(b - a) / denom:10.2e}")


if __name__ == "__main__":
    _demo()
