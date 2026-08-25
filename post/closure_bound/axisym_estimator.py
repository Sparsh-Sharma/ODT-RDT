"""
LOS -> transverse-plane estimator for the rapid pressure-strain closure.

Implements the machinery of notes/los_estimator_note.tex in the conventions of
the JFM-2026-1781 manuscript, Section 5:

  line along e_2;  kperp = sqrt(k1^2 + k3^2);  k^2 = k2^2 + kperp^2;
  angular variable x = k2^2 / k^2;
  axisymmetric solenoidal spectrum tensor (assumption A1)
      Phi_mn = a(k2,kperp) P_mn + c(k2,kperp) xi_m xi_n ,
      P_mn = delta_mn - k_m k_n / k^2 ,  xi = e2 - mu * k/|k| ,  mu^2 = x,
  realizability cone:  a >= 0  and  a + (1-x) c >= 0.

Line-of-sight (LOS) data: the one-dimensional spectra
      phi_mn(k2) = int int Phi_mn dk1 dk3 .
Under A1 the forward map is (note eqs. (5)-(7)):
      phi_11 = phi_33 = 2*pi * int [ a*(1+x)/2 + c*x*(1-x)/2 ] kperp dkperp
      phi_22          = 2*pi * int [ a*(1-x)   + c*(1-x)^2   ] kperp dkperp
      phi_mn (m != n) = 0.

Target functional: the collapsed rapid pressure-strain term for plane strain
A = diag(1/2, -1/2, 0), manuscript eqs. (5.7)-(5.10), with kernels
      g1 = a*(1/8 + 3/4 x - 7/8 x^2) + c * 7/8 * x (1-x)^2
      g2 = -3/2 x * [ a (1-x) + c (1-x)^2 ]
      g3 = a*(-1/8 + 3/4 x - 5/8 x^2) + c * 5/8 * x (1-x)^2 .

Normalization convention (pinned by the exact isotropic limit): we define
      Pi_nn = 2 * int_{R^3} g_n(x)[a, c] d^3kappa
            = 4 * int_0^inf dk2 int_0^inf 2*pi*kperp dkperp g_n(x)[a, c]
(the second equality uses k2 -> -k2 symmetry). With a_iso = E(k)/(4*pi*k^2)
and c = 0 this yields exactly Pi = (4/5) k_t S_ij, i.e. (2/5, -2/5, 0) k_t
for the plane strain above -- the Crow (1968) result the manuscript quotes.

Estimator: at each k2 the data constraints, the realizability cone, and the
objective are all linear in (a, c)(kperp), so sharp global bounds on
pi_n(k2) = 4 * 2*pi * int g_n kperp dkperp follow from two small linear
programs (scipy linprog, HiGHS). Everything decomposes per k2.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import linprog

TWO_PI = 2.0 * np.pi
# Pi_nn = PI_PREFACTOR * int_0^inf dk2 [ 2*pi int_0^inf g_n kperp dkperp ]
# (factor 2 for the k2<0 half-space, factor 2 from the Poisson-source
#  convention of the manuscript; pinned by test_exact_pi_isotropic).
PI_PREFACTOR = 4.0


# ----------------------------------------------------------------------------
# kernels (manuscript eqs. 5.8-5.10), coefficients of the two spectral scalars
# ----------------------------------------------------------------------------

def g_kernels_a(x):
    """Coefficients of the isotropic-like scalar a in (g1, g2, g3)."""
    x = np.asarray(x, dtype=float)
    g1 = 0.125 + 0.75 * x - 0.875 * x**2
    g2 = -1.5 * x * (1.0 - x)
    g3 = -0.125 + 0.75 * x - 0.625 * x**2
    return g1, g2, g3


def g_kernels_c(x):
    """Coefficients of the anisotropy scalar c in (g1, g2, g3)."""
    x = np.asarray(x, dtype=float)
    g1 = 0.875 * x * (1.0 - x) ** 2
    g2 = -1.5 * x * (1.0 - x) ** 2
    g3 = 0.625 * x * (1.0 - x) ** 2
    return g1, g2, g3


# ----------------------------------------------------------------------------
# grids and quadrature
# ----------------------------------------------------------------------------

@dataclass
class Grid:
    """Log-spaced quadrature grids for the (k2 > 0, kperp > 0) quarter-plane."""

    k2: np.ndarray      # (N,)  positive half-axis
    kperp: np.ndarray   # (M,)
    w2: np.ndarray      # (N,)  trapezoid weights for the k2 integral
    wperp: np.ndarray   # (M,)  trapezoid weights INCLUDING the 2*pi*kperp factor

    @classmethod
    def make(cls, k_min=1e-3, k_max=1e3, n_k2=160, n_kperp=240):
        k2 = np.geomspace(k_min, k_max, n_k2)
        kperp = np.geomspace(k_min, k_max, n_kperp)
        return cls(k2=k2, kperp=kperp,
                   w2=_trapz_weights(k2),
                   wperp=TWO_PI * kperp * _trapz_weights(kperp))

    def x(self, k2_val):
        """Angular variable x = k2^2/(k2^2 + kperp^2) along the kperp grid."""
        return k2_val**2 / (k2_val**2 + self.kperp**2)


def _trapz_weights(z):
    w = np.zeros_like(z)
    w[1:-1] = 0.5 * (z[2:] - z[:-2])
    w[0] = 0.5 * (z[1] - z[0])
    w[-1] = 0.5 * (z[-1] - z[-2])
    return w


# ----------------------------------------------------------------------------
# model spectra (for tests / synthetic data)
# ----------------------------------------------------------------------------

def vk_energy_spectrum(k, L=1.0, amplitude=1.0):
    """von Karman energy spectrum shape E(k) ~ (kL)^4 / (1+(kL)^2)^(17/6)."""
    kl = k * L
    return amplitude * L * kl**4 / (1.0 + kl**2) ** (17.0 / 6.0)


def a_isotropic(k2_val, kperp, L=1.0, amplitude=1.0):
    """Isotropic scalar a(k2,kperp) = E(k)/(4 pi k^2) on the kperp grid."""
    k = np.sqrt(k2_val**2 + kperp**2)
    return vk_energy_spectrum(k, L, amplitude) / (4.0 * np.pi * k**2)


# ----------------------------------------------------------------------------
# forward map and exact functionals (given a, c on the grid)
# ----------------------------------------------------------------------------

def line_spectra(grid: Grid, a_fn, c_fn=None):
    """phi_11 (= phi_33) and phi_22 on grid.k2 from scalar fields a, c.

    a_fn, c_fn: callables (k2_val, kperp_array) -> array. c_fn=None means c=0.
    """
    n = grid.k2.size
    phi11 = np.empty(n)
    phi22 = np.empty(n)
    for i, k2v in enumerate(grid.k2):
        x = grid.x(k2v)
        a = a_fn(k2v, grid.kperp)
        c = c_fn(k2v, grid.kperp) if c_fn is not None else np.zeros_like(a)
        phi11[i] = np.sum(grid.wperp * (a * 0.5 * (1.0 + x)
                                        + c * 0.5 * x * (1.0 - x)))
        phi22[i] = np.sum(grid.wperp * (a * (1.0 - x)
                                        + c * (1.0 - x) ** 2))
    return phi11, phi22


def kinetic_energy(grid: Grid, a_fn, c_fn=None):
    """k_t = (1/2) int Phi_mm d^3kappa, full space = 2 x (k2>0 half-space)."""
    kt = 0.0
    for i, k2v in enumerate(grid.k2):
        x = grid.x(k2v)
        a = a_fn(k2v, grid.kperp)
        c = c_fn(k2v, grid.kperp) if c_fn is not None else np.zeros_like(a)
        trace = 2.0 * a + c * (1.0 - x)
        kt += grid.w2[i] * np.sum(grid.wperp * trace)
    return 0.5 * 2.0 * kt


def pi_exact(grid: Grid, a_fn, c_fn=None):
    """Exact (Pi_11, Pi_22, Pi_33) for known scalars a, c on the grid."""
    pi = np.zeros(3)
    for i, k2v in enumerate(grid.k2):
        x = grid.x(k2v)
        a = a_fn(k2v, grid.kperp)
        c = c_fn(k2v, grid.kperp) if c_fn is not None else np.zeros_like(a)
        ga = g_kernels_a(x)
        gc = g_kernels_c(x)
        for n in range(3):
            pi[n] += grid.w2[i] * np.sum(grid.wperp * (ga[n] * a + gc[n] * c))
    return PI_PREFACTOR * pi


# ----------------------------------------------------------------------------
# the LP estimator
# ----------------------------------------------------------------------------

@dataclass
class BoundResult:
    """Per-k2 sharp bounds on the rapid pressure-strain density pi_n(k2)."""

    k2: np.ndarray            # (N,)
    lower: np.ndarray         # (3, N)  min of pi_n density at each k2
    upper: np.ndarray         # (3, N)  max
    status_ok: np.ndarray     # (N,) bool, both LPs solved at this k2

    def integrated(self, w2):
        """Integrated bounds (Pi^-, Pi^+), each shape (3,)."""
        lo = (self.lower * w2).sum(axis=1)
        hi = (self.upper * w2).sum(axis=1)
        return lo, hi


def _slope_rows(m, offset, s, n_var, wperp):
    """Rows enforcing the log-Lipschitz cap |d ln z / d ln kperp| <= lambda
    on the PHYSICAL field z_phys = const * z_var / wperp, for a nonnegative
    variable block starting at `offset` (stage b).

    Physical form: (1-s) z_{j+1} <= (1+s) z_j and mirrored, with
    s_j = tanh(lambda dln kperp_j / 2), which is the exact discrete cap
    z_{j+1}/z_j <= exp(lambda dln kperp). In scaled variables each row is
    multiplied by w_j w_{j+1} / (w_j + w_{j+1}) so all coefficients are O(1).
    """
    rows = np.zeros((2 * (m - 1), n_var))
    j = np.arange(m - 1)
    p = wperp[:-1] / (wperp[:-1] + wperp[1:])
    q = wperp[1:] / (wperp[:-1] + wperp[1:])
    # growth cap: (1-s) p_j z_{j+1} - (1+s) q_j z_j <= 0
    rows[2 * j, offset + j] = -(1.0 + s) * q
    rows[2 * j, offset + j + 1] = (1.0 - s) * p
    # decay cap: (1-s) q_j z_j - (1+s) p_j z_{j+1} <= 0
    rows[2 * j + 1, offset + j] = (1.0 - s) * q
    rows[2 * j + 1, offset + j + 1] = -(1.0 + s) * p
    return rows


def bound_pi_at_k2(k2_val, kperp, wperp, phi11, phi22, component,
                   data_band=0.0, slope_cap=None, polarization_cap=None):
    """Sharp [min, max] of the pi_n density at one k2 from line data alone.

    Density convention: pi_n(k2) = PI_PREFACTOR * 2*pi*int g_n kperp dkperp,
    so that int pi_n dk2 over the k2>0 grid equals Pi_nn.

    Variables z = [a_1..a_M, cp_1..cp_M, cm_1..cm_M] on the kperp grid, all
    nonnegative, with c = cp - cm (split so the sign-indefinite anisotropy
    scalar admits the relative-slope constraint below).
    Constraints:
      (D1) sum wperp * (a(1+x)/2 + c x(1-x)/2) = phi11   (+/- band)
      (D2) sum wperp * (a(1-x)   + c (1-x)^2 ) = phi22   (+/- band)
      (R)  a_j >= 0,  a_j + (1-x_j) c_j >= 0.
      (S)  stage b, if slope_cap = lambda is given: log-Lipschitz cap
           |d ln z / d ln kperp| <= lambda on each of a, cp, cm, imposed as
           |z_{j+1}-z_j| <= s_j (z_j + z_{j+1}), s_j = tanh(lambda dlnk_j / 2)
           (exact discrete form of the cap for geometric growth). This is
           scale-free, keeps the LP linear, and excludes the bang-bang
           extremals of the raw bound. The von Karman truth needs
           lambda >= 11/3 (tail a ~ kperp^{-11/3}).
    data_band: relative half-width of an inequality band replacing the
    equalities (0 -> exact equalities; use >0 for noisy/measured data).
    """
    m = kperp.size
    n_var = 3 * m
    x = k2_val**2 / (k2_val**2 + kperp**2)

    # scaled variables u_j = wperp_j * (field)_j / phi22: every constraint
    # block then has O(1) coefficients and O(1) rhs (HiGHS fails outright on
    # the unscaled problem once slope rows couple columns spanning decades);
    # the physical objective is recovered by the factor PI_PREFACTOR * phi22.
    amp = max(abs(phi22), 1e-300)
    ratio11 = phi11 / amp

    ga = g_kernels_a(x)[component]
    gc = g_kernels_c(x)[component]
    obj = np.concatenate([ga, gc, -gc])
    obj_back = PI_PREFACTOR * amp

    def data_row(wa, wc):
        return np.concatenate([wa, wc, -wc])

    row11 = data_row(0.5 * (1.0 + x), 0.5 * x * (1.0 - x))
    row22 = data_row(1.0 - x, (1.0 - x) ** 2)

    # realizability: -(u^a_j + (1-x_j)(u^cp_j - u^cm_j)) <= 0
    a_real = np.zeros((m, n_var))
    j = np.arange(m)
    a_real[j, j] = -1.0
    a_real[j, m + j] = -(1.0 - x)
    a_real[j, 2 * m + j] = (1.0 - x)
    b_real = np.zeros(m)

    ub_blocks, ub_rhs = [a_real], [b_real]

    if slope_cap is not None:
        s = np.tanh(0.5 * slope_cap * np.diff(np.log(kperp)))
        for block in range(3):
            ub_blocks.append(_slope_rows(m, block * m, s, n_var, wperp))
            ub_rhs.append(np.zeros(2 * (m - 1)))

    if polarization_cap is not None:
        # |c_j| <= gamma a_j pointwise (gamma=0: isotropic polarization).
        # In scaled variables the common wperp/phi22 factor cancels.
        gam = float(polarization_cap)
        pol = np.zeros((2 * m, n_var))
        pol[j, j] = -gam
        pol[j, m + j] = 1.0
        pol[j, 2 * m + j] = -1.0
        pol[m + j, j] = -gam
        pol[m + j, m + j] = -1.0
        pol[m + j, 2 * m + j] = 1.0
        ub_blocks.append(pol)
        ub_rhs.append(np.zeros(2 * m))

    if data_band > 0.0:
        band11 = data_band * abs(ratio11)
        band22 = data_band
        ub_blocks += [row11[None, :], -row11[None, :],
                      row22[None, :], -row22[None, :]]
        ub_rhs += [[ratio11 + band11], [-(ratio11 - band11)],
                   [1.0 + band22], [-(1.0 - band22)]]
        a_eq, b_eq = None, None
    else:
        a_eq = np.vstack([row11, row22])
        b_eq = np.array([ratio11, 1.0])

    a_ub = np.vstack(ub_blocks)
    b_ub = np.concatenate([np.atleast_1d(r) for r in ub_rhs])
    bounds = [(0.0, None)] * n_var

    out = {}
    for sign, key in ((+1.0, "min"), (-1.0, "max")):
        res = linprog(sign * obj, A_ub=a_ub, b_ub=b_ub,
                      A_eq=a_eq, b_eq=b_eq, bounds=bounds, method="highs")
        out[key] = sign * obj_back * res.fun if res.status == 0 else np.nan
        out[key + "_ok"] = res.status == 0
    return out["min"], out["max"], out["min_ok"] and out["max_ok"]


def bound_pi(grid: Grid, phi11, phi22, data_band=0.0, slope_cap=None,
             polarization_cap=None):
    """Run the per-k2 LPs for all three components over the whole k2 grid."""
    n = grid.k2.size
    lower = np.full((3, n), np.nan)
    upper = np.full((3, n), np.nan)
    ok = np.zeros(n, dtype=bool)
    for i, k2v in enumerate(grid.k2):
        all_ok = True
        for comp in range(3):
            lo, hi, solved = bound_pi_at_k2(k2v, grid.kperp, grid.wperp,
                                            phi11[i], phi22[i], comp,
                                            data_band=data_band,
                                            slope_cap=slope_cap,
                                            polarization_cap=polarization_cap)
            lower[comp, i], upper[comp, i] = lo, hi
            all_ok &= solved
        ok[i] = all_ok
    return BoundResult(k2=grid.k2, lower=lower, upper=upper, status_ok=ok)


# ----------------------------------------------------------------------------
# isotropic-limit calibration (the first real number)
# ----------------------------------------------------------------------------

def isotropic_calibration(L=1.0, amplitude=1.0, slope_cap=None,
                          polarization_cap=None, **grid_kw):
    """Bounds vs exact Pi for an isotropic von Karman spectrum.

    Returns dict with k_t, exact Pi, integrated bounds, and the normalized
    indeterminacy I_n = (Pi+ - Pi-)/(2 k_t S) with S = 1/2 (plane strain).
    slope_cap: stage-b log-Lipschitz cap lambda (None = no cap).
    polarization_cap: stage-b cap gamma in |c| <= gamma a (None = no cap).
    """
    grid = Grid.make(**grid_kw)
    a_fn = lambda k2v, kp: a_isotropic(k2v, kp, L, amplitude)
    phi11, phi22 = line_spectra(grid, a_fn)
    kt = kinetic_energy(grid, a_fn)
    exact = pi_exact(grid, a_fn)
    bounds = bound_pi(grid, phi11, phi22, slope_cap=slope_cap,
                      polarization_cap=polarization_cap)
    lo, hi = bounds.integrated(grid.w2)
    indeterminacy = (hi - lo) / (2.0 * kt * 0.5)
    return {"k_t": kt, "pi_exact": exact, "pi_lower": lo, "pi_upper": hi,
            "indeterminacy": indeterminacy, "bounds": bounds, "grid": grid}


def stage_b_sweep(cases, L=1.0, amplitude=1.0, **grid_kw):
    """Stage-b sweep over constraint combinations.

    cases: list of (label, dict) where the dict holds slope_cap and/or
    polarization_cap. Returns list of result dicts with 'bracketing' =
    does [Pi-, Pi+] still contain the exact isotropic Pi.
    """
    rows = []
    for label, kw in cases:
        cal = isotropic_calibration(L=L, amplitude=amplitude, **kw, **grid_kw)
        tol = 5e-3 * cal["k_t"]
        bracketing = bool(
            np.all(cal["pi_lower"] <= cal["pi_exact"] + tol)
            and np.all(cal["pi_exact"] - tol <= cal["pi_upper"]))
        rows.append({"label": label, "bracketing": bracketing, **cal})
    return rows


def _print_calibration(cal, kt):
    print("component   exact/k_t    lower/k_t    upper/k_t    I_n")
    for n, name in enumerate(("11", "22", "33")):
        print(f"Pi_{name}      {cal['pi_exact'][n]/kt:+.4f}      "
              f"{cal['pi_lower'][n]/kt:+.4f}      "
              f"{cal['pi_upper'][n]/kt:+.4f}      "
              f"{cal['indeterminacy'][n]:.4f}")


if __name__ == "__main__":
    cases = [
        ("raw (stage a)", {}),
        ("slope lambda=4", dict(slope_cap=4.0)),
        ("polarization gamma=1", dict(polarization_cap=1.0)),
        ("polarization gamma=0.3", dict(polarization_cap=0.3)),
        ("polarization gamma=0", dict(polarization_cap=0.0)),
        ("gamma=0 + lambda=4", dict(polarization_cap=0.0, slope_cap=4.0)),
    ]
    rows = stage_b_sweep(cases, n_k2=48, n_kperp=160, k_min=1e-2, k_max=1e2)
    kt = rows[0]["k_t"]
    print(f"k_t = {kt:.6f}   (isotropic exact Pi/k_t = +2/5, -2/5, 0)\n")
    for row in rows:
        print(f"--- {row['label']}   bracketing: {row['bracketing']}")
        _print_calibration(row, kt)
        print()
