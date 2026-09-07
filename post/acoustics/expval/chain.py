#!/usr/bin/env python3
"""Four-way absolute Amiet prediction chain for the experimental
validation cases (Part III extension): classical Amiet + von Karman,
Amiet + Liepmann, Amiet + standard (strain-coupled baseline) ODT, and
Amiet + fixed ODT (size-capped relaxation clock, RCS1).

Acoustic stage: the compact lift-dipole far-field PSD with the exact
Sears response - the normalisation-audited len_vsdb form (verbatim from
amiet_reference.far_field_psd_dipole) - with the spectrum-consistent
spanwise correlation length l_y = pi*Phi(Kx,0)/Int Phi dKy for EVERY
inflow model, so the four predictions differ only in the inflow.

Inflow stage:
  vK / Liepmann : analytic isotropic 2D upwash spectra at the measured
                  (u'_w, Lambda_f).
  ODT (std/fix) : the exact-RDT-distorted-vK representation fitted to
                  the model's strained line ensembles (A0, ke per
                  strain; rdt_kernel), bridged to physical units by
                  kappa_e = 0.7468/Lambda_f and the e=0 upwash
                  variance = u'_w^2; level change from the measured
                  ensemble variance ratio; evaluated at the case's
                  effective distortion strain e_eff (potential-flow
                  stagnation-line strain of the section, truncated at
                  the eddy scale - see e_eff_ellipse).

Response-side thickness (NACA sections only): len_vsdb ThickAirfoil
Gershfeld surrogate G_thick = (1+(k r_LE/b)^2)^(-5/6), applied
identically to all four chains; quantitative metrics are restricted to
f < U/t_A where G_thick is a small correction.
"""
import os
import re
import sys

import numpy as np
from scipy.special import fresnel, jv, yv

HERE = os.path.dirname(os.path.abspath(__file__))
LENOISE = os.path.normpath(os.path.join(
    HERE, "..", "..", "closure_bound", "lenoise"))
sys.path.insert(0, LENOISE)
import rdt_kernel as rk  # noqa: E402

RHO0, C0, PREF = 1.225, 340.0, 20e-6
CVK = 0.746834

# ----------------------------------------------------------------- acoustics


def sears(k):
    """Exact Sears function (len_vsdb gust_response)."""
    k = np.asarray(k, float)
    H1 = jv(1, k) - 1j * yv(1, k)
    H0 = jv(0, k) - 1j * yv(0, k)
    C = H1 / (H1 + 1j * H0)
    return C * (jv(0, k) - 1j * jv(1, k)) + 1j * jv(1, k)


def far_field_psd_dipole(x, omega, S, phi1d, chord, U, span, M, ell_y):
    """len_vsdb amiet_reference.far_field_psd_dipole (verbatim);
    phi1d = 1-D upwash wavenumber spectrum (two-sided, int dKx = w'^2)."""
    x1, x2, x3 = x
    sigma2 = x1 ** 2 + (1 - M ** 2) * (x2 ** 2 + x3 ** 2)
    dipole2 = (omega * x3 / (4 * np.pi * C0 * sigma2)) ** 2
    lift_per_w2 = (np.pi * RHO0 * U * chord) ** 2 * np.abs(S) ** 2
    return dipole2 * lift_per_w2 * (phi1d / U) * span * ell_y


def _E(xi):
    """Fresnel combination E(xi) = int_0^xi e^{it}/sqrt(2 pi t) dt
    = C(sqrt(2 xi/pi)) + i S(sqrt(2 xi/pi))."""
    z = np.sqrt(2.0 * np.asarray(xi, float) / np.pi)
    S, C = fresnel(z)
    return C + 1j * S


def amiet_L(x1_S0, k1s, mu, M, beta, k2s=0.0):
    """Amiet's chordwise non-compact LE transfer function L = L1 + L2
    (supercritical gusts, Bampanis-Roger-Moreau 2022 eqs (3)-(4));
    k1s = omega c/(2U) (dimensionless chordwise wavenumber), mu =
    k1s M/beta^2, x1_S0 = x1/S0.  Valid mu > ~0.4."""
    k1s = np.asarray(k1s, float)
    mu = np.asarray(mu, float)
    kap = np.sqrt(np.maximum(mu ** 2 - (k2s / beta) ** 2, 1e-12))
    th2 = mu * (M - x1_S0) - np.pi / 4
    th3 = kap + mu * x1_S0
    th4 = kap - mu * x1_S0
    pre = k1s + beta ** 2 * kap
    L1 = -(1.0 / np.pi) * np.sqrt(2.0 / (pre * th4)) \
        * np.exp(-1j * th2) * _E(2 * th4)
    L2 = np.exp(-1j * th2) / (np.pi * np.sqrt(2 * np.pi * pre * th4)) \
        * (1j * (1 - np.exp(2j * th4))
           - (1 + 1j) * (_E(4 * kap)
                         - np.exp(2j * th4) * np.sqrt(2 * kap / th3)
                         * _E(2 * th3)))
    return L1 + L2


def g_thick(kred, t_over_c):
    """Gershfeld thickness admittance (len_vsdb ThickAirfoil.admittance),
    kred = omega*b/U, r_LE/b = 2*1.1019*(t/c)^2."""
    rle_b = 2.0 * 1.1019 * t_over_c ** 2
    return (1.0 + (np.asarray(kred, float) * rle_b) ** 2) ** (-5.0 / 6.0)


# ------------------------------------------------------------ inflow spectra


class Spectrum2D:
    """Interface: phi(Kx, Ky) 2D upwash spectrum [m^4/s^2] (two-sided)."""

    def phi(self, kx, ky):
        raise NotImplementedError

    def phi1d(self, kx, nky=160, half=9.0):
        """Int phi dKy and l_y = pi*phi(Kx,0)/Int phi dKy, per Kx."""
        kx = np.atleast_1d(np.asarray(kx, float))
        ky = self._kygrid(nky, half)
        P = np.array([self.phi(k, ky) for k in kx])
        I = np.trapezoid(P, ky, axis=1)
        p0 = np.array([float(self.phi(k, np.array([0.0]))[0]) for k in kx])
        return I, np.pi * p0 / I

    def _kygrid(self, n, half):
        x = self.kscale * np.exp(np.linspace(-half, np.log(30.0), n // 2))
        return np.concatenate([-x[::-1], x])


class VonKarman2D(Spectrum2D):
    """Isotropic vK 2D upwash spectrum (Amiet 1975):
    Phi_ww = (4/(9 pi)) (w2/ke^2) (kx^2+ky^2)/ke^2 / (1+(kx^2+ky^2)/ke^2)^(7/3)
    normalised numerically to integrate to w2 exactly."""

    def __init__(self, w2, Lam):
        self.ke = CVK / Lam
        self.kscale = self.ke
        self.w2 = w2
        self._A = 1.0
        self._A = w2 / self._variance()

    def phi(self, kx, ky):
        kh2 = (np.asarray(kx, float) ** 2 + np.asarray(ky, float) ** 2) \
            / self.ke ** 2
        return self._A * (4.0 / (9 * np.pi)) / self.ke ** 2 \
            * kh2 / (1 + kh2) ** (7.0 / 3.0)

    def _variance(self, n=200, half=9.0):
        g = self._kygrid(n, half)
        rows = np.array([np.trapezoid(self.phi(kx, g), g) for kx in g])
        return np.trapezoid(rows, g)


class Liepmann2D(Spectrum2D):
    """Isotropic Liepmann 2D upwash spectrum:
    Phi_ww = (3 w2 L^2/(4 pi)) L^2(kx^2+ky^2) / (1+L^2(kx^2+ky^2))^(5/2),
    normalised numerically (guards the analytic constant)."""

    def __init__(self, w2, Lam):
        self.L = Lam
        self.kscale = 1.0 / Lam
        self.w2 = w2
        self._A = 1.0
        self._A = w2 / self._variance()

    def phi(self, kx, ky):
        s = self.L ** 2 * (np.asarray(kx, float) ** 2
                           + np.asarray(ky, float) ** 2)
        return self._A * (3 * self.w2 * self.L ** 2 / (4 * np.pi)) \
            * s / (1 + s) ** 2.5 / self.w2

    def _variance(self, n=200, half=9.0):
        g = self._kygrid(n, half)
        rows = np.array([np.trapezoid(self.phi(kx, g), g) for kx in g])
        return np.trapezoid(rows, g)


def read_fits(tag):
    fits = {}
    for line in open(os.path.join(LENOISE,
                                  f"rdt_family_fit_table_{tag}.txt")):
        m = re.match(r"\s*([\d.]+)\s+([\d.]+)\s+([\d.eE+-]+)\s+[\d.]+",
                     line)
        if m:
            fits[float(m.group(1))] = (float(m.group(2)),
                                       float(m.group(3)))
    return fits


class ODTRep2D(Spectrum2D):
    """The RDT-distorted-vK representation of an ODT strained ensemble,
    in PHYSICAL units, at arbitrary strain e (parameters interpolated
    log-linearly between the fitted strains); amplitude bridged so the
    e=0 member carries variance w2 at length scale Lam, and the level
    at e comes from the measured ensemble variance ratio r22(e)."""

    def __init__(self, w2, Lam, e, fit_table, spectra_npz):
        self.e = float(e)
        fits = fit_table
        es = np.array(sorted(fits))
        kes = np.array([fits[x][0] for x in es])
        ke_m = float(np.exp(np.interp(e, es, np.log(kes))))
        ke0_m = fits[0.0][0]
        self.Lunit = ke0_m / (CVK / Lam)      # metres per model unit
        self.ke_m = ke_m
        self.kscale = ke_m / self.Lunit
        # measured level ratio (dumps ordered e=0,0.5,1,1.5,2)
        sp = np.load(spectra_npz, allow_pickle=True)

        def r22(ee):
            i = int(round(ee / 0.5))
            return np.trapezoid(sp[f"med_E2_e{i}"], sp[f"k2_e{i}"])

        eg = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
        rg = np.array([r22(x) for x in eg])
        lev = float(np.exp(np.interp(e, eg, np.log(rg))) / rg[0])
        # normalise the representation to unit variance, then w2*lev
        self._A = 1.0
        self._A = w2 * lev / self._variance()

    def phi(self, kx, ky):
        """model kernel evaluated at physical (kx, ky)."""
        kxm = np.asarray(kx, float) * self.Lunit
        kym = np.asarray(ky, float) * self.Lunit
        k2 = self.ke_m * np.exp(np.linspace(-9.0, np.log(30.0), 120))
        k2 = np.concatenate([-k2[::-1], k2])
        out = np.empty(np.broadcast(kxm, kym).shape, float)
        kxb, kyb = np.broadcast_arrays(kxm, kym)
        flat = out.reshape(-1)
        for i, (a, b) in enumerate(zip(kxb.reshape(-1), kyb.reshape(-1))):
            P22, _ = rk.rdt_phi_components(a / self.ke_m, k2 / self.ke_m,
                                           b / self.ke_m, self.e)
            flat[i] = np.trapezoid(P22, k2)
        return out * self._A

    def _variance(self, n=90, half=9.0):
        g = self._kygrid(n, half)
        rows = np.array([np.trapezoid(self.phi(kx, g), g) for kx in g])
        return np.trapezoid(rows, g)


# ------------------------------------------------------- effective distortion


def e_eff_ellipse(chord, thick, xstar):
    """Accumulated plane strain along the stagnation streamline of a
    potential-flow ellipse (semi-axes a=chord/2, b=thick/2), truncated
    where the eddy centre stands xstar ahead of the leading edge:
    e_eff = ln(U_inf / u(xstar)).  Joukowski map, evaluated exactly."""
    a, b = chord / 2.0, thick / 2.0
    r = (a + b) / 2.0                          # mapped cylinder radius
    m = (a ** 2 - b ** 2) / 4.0                # z = zeta + m/zeta
    x = a + xstar                              # physical point on axis
    # invert z = zeta + m/zeta on the real axis (zeta >= r)
    zeta = 0.5 * (x + np.sqrt(x ** 2 - 4 * m))
    w = (1.0 - r ** 2 / zeta ** 2) / (1.0 - m / zeta ** 2)
    return float(np.log(1.0 / max(w, 1e-9)))


# ---------------------------------------------------------------- prediction


def predict_spl(f, U, w2, Lam, chord, span, obs, model, t_over_c=None,
                fit_tag=None, e_eff=0.0):
    """SPL(f) dB/Hz re 20 uPa for one inflow model.
    model in {'vk','liepmann','odt_std','odt_fix'}; obs = (x1,x2,x3);
    fit_tag e.g. 'S20' (standard) - the fixed model reads the *_RCS1
    tables when model='odt_fix'."""
    omega = 2 * np.pi * np.asarray(f, float)
    Kx = omega / U
    M = U / C0
    kred = omega * (chord / 2) / U
    S = sears(kred)
    if model == "vk":
        sp = VonKarman2D(w2, Lam)
    elif model == "liepmann":
        sp = Liepmann2D(w2, Lam)
    elif model == "odt_std":
        sp = ODTRep2D(w2, Lam, e_eff, read_fits(fit_tag),
                      os.path.join(LENOISE, f"spectra_gateA_{fit_tag}.npz"))
    elif model == "odt_fix":
        sp = ODTRep2D(w2, Lam, e_eff, read_fits(fit_tag + "_RCS1"),
                      os.path.join(LENOISE,
                                   f"spectra_gateA_{fit_tag}_RCS1.npz"))
    else:
        raise ValueError(model)
    # Amiet's large-aspect-ratio far-field PSD (Bampanis et al. 2022
    # eq. (2)): chordwise non-compact response, spanwise wavenumber
    # selected by the observer direction (k2 = k x2/S0).
    x1, x2, x3 = obs
    beta = np.sqrt(1 - M ** 2)
    S0 = np.sqrt(x1 ** 2 + beta ** 2 * (x2 ** 2 + x3 ** 2))
    kac = omega / C0
    k2sel = kac * x2 / S0
    mu = kred * M / beta ** 2
    L = amiet_L(x1 / S0, kred, mu, M, beta, k2s=k2sel * chord / 2)
    phi2d = np.array([float(sp.phi(kx, np.array([kk]))[0])
                      for kx, kk in zip(Kx, np.broadcast_to(
                          k2sel, Kx.shape))])
    L2 = np.abs(L) ** 2
    if t_over_c is not None:
        L2 = L2 * g_thick(kred, t_over_c)
    Spp = (kac * RHO0 * chord * x3 / (2 * S0 ** 2)) ** 2 * np.pi * U \
        * (span / 2.0) * phi2d * L2
    # Spp is two-sided in omega (Amiet's convention): one-sided per-Hz
    # level needs Gpp = 2*Spp and the 2*pi Hz<->rad/s factor (CR-2733
    # p.22 states exactly this bookkeeping).
    return 10 * np.log10(Spp * 4 * np.pi / PREF ** 2)
