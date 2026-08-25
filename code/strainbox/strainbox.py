"""strainbox: serial pseudo-spectral DNS of homogeneous turbulence under
uniform irrotational mean strain, in the Rogallo (1981) deforming frame.

Scope: the physics prototype for cases/dns_strainedBox/SCOPING.md — correct
strain coupling, validated against an independent per-mode rapid-distortion
integrator (test_strainbox.py). Production MPI port comes after this is
proven.

Formulation
-----------
Mean flow U_i = A_ij x_j with A = S diag(a1, a2, a3), sum a_i = 0
(irrotational plane strain: a = (1/2, -1/2, 0)). In the deforming frame the
physical wavevector of a mode with label k0 is

    k_i(t) = k0_i * exp(-a_i * e),   e = S t   (accumulated strain),

and the fluctuation spectrum obeys (incompressible, solenoidal)

    d u^_i/dt = -A_ij u^_j + 2 (k_i k_m / k^2) A_mn u^_n      [RDT part]
                - P_il(k) N^_l(u)                             [nonlinear]
                - nu k^2 u^_i                                 [viscous]

with P_il = delta_il - k_i k_l / k^2 and N^ the Fourier transform of
u_j du_i/dx_j evaluated pseudo-spectrally with the CURRENT k(t) and 2/3
dealiasing. Dropping the nonlinear and viscous terms gives exact linear RDT
— the acceptance test compares the full solver at large S k/eps against an
independent integration of that linear system.

Conventions: box (2 pi)^3 at t=0, integer k0, numpy rfftn layout
(kx, ky half-axis last? -> we use axes (0,1,2) = (x1, x2, x3), rfft over
axis 2). Velocities real on the N^3 grid; spectra normalized so that
sum_k |u^|^2 (with rfft weights) = <u_i u_i>/... see energy().
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class StrainBox:
    n: int = 64
    nu: float = 0.0
    smag: float = 0.0                      # strain rate S
    a_dir: tuple = (0.5, -0.5, 0.0)        # A = S diag(a_dir)
    dealias: bool = True
    seed: int = 12345
    linear_only: bool = False              # drop nonlinear term (RDT mode)

    def __post_init__(self):
        n = self.n
        k1 = np.fft.fftfreq(n, 1.0 / n)            # integer wavenumbers
        self.k0 = np.array(np.meshgrid(k1, k1, k1[:n // 2 + 1],
                                       indexing="ij"))     # (3,n,n,n/2+1)
        self.e = 0.0
        self.t = 0.0
        # 2/3-rule mask on the LABEL grid (k0): modes are carried, not
        # remeshed, so the mask is fixed
        kmax = n // 3
        self.mask = ((np.abs(self.k0[0]) <= kmax)
                     & (np.abs(self.k0[1]) <= kmax)
                     & (np.abs(self.k0[2]) <= kmax))
        # rfft double-count weights for spectral sums (last axis)
        w = np.full(self.k0.shape[1:], 2.0)
        w[..., 0] = 1.0
        if n % 2 == 0:
            w[..., -1] = 1.0
        self.specw = w
        self.uh = np.zeros((3,) + self.k0.shape[1:], dtype=np.complex128)

    # ------------------------------------------------------------------ #
    def k_phys(self, e=None):
        """Physical wavevector at accumulated strain e (default: current)."""
        e = self.e if e is None else e
        stretch = np.exp(-np.asarray(self.a_dir) * e)
        return self.k0 * stretch[:, None, None, None]

    def project(self, uh, k):
        k2 = (k ** 2).sum(axis=0)
        k2[0, 0, 0] = 1.0
        div = (k * uh).sum(axis=0) / k2
        return uh - k * div[None, ...]

    # ------------------------------------------------------------------ #
    def init_isotropic(self, spectrum, kt_target=None):
        """Random-phase solenoidal field with E(k) = spectrum(k) at t=0."""
        rng = np.random.default_rng(self.seed)
        shape = self.k0.shape[1:]
        uh = (rng.standard_normal((3,) + shape)
              + 1j * rng.standard_normal((3,) + shape))
        k = self.k0
        kmag = np.sqrt((k ** 2).sum(axis=0))
        kmag[0, 0, 0] = 1.0
        uh = self.project(uh, k)
        # scale shells to the target spectrum
        with np.errstate(divide="ignore", invalid="ignore"):
            target = spectrum(kmag)
        target[kmag < 0.5] = 0.0
        cur2 = (np.abs(uh) ** 2).sum(axis=0) * self.specw
        shell = np.rint(kmag).astype(int)
        nsh = shell.max() + 1
        e_cur = np.bincount(shell.ravel(), cur2.ravel(), minlength=nsh)
        e_tgt = np.zeros(nsh)
        for s in range(1, nsh):
            e_tgt[s] = spectrum(float(s))
        scale = np.zeros(nsh)
        nz = e_cur > 0
        scale[nz] = np.sqrt(e_tgt[nz] / e_cur[nz])
        uh *= scale[shell][None, ...]
        uh[:, 0, 0, 0] = 0.0
        uh *= self.mask[None, ...]
        # enforce Hermitian symmetry via round-trip
        for c in range(3):
            uh[c] = np.fft.rfftn(np.fft.irfftn(uh[c], s=(self.n,) * 3, axes=(0, 1, 2)))
        self.uh = self.project(uh, k)
        if kt_target is not None:
            self.uh *= np.sqrt(kt_target / self.kinetic_energy())

    # ------------------------------------------------------------------ #
    def kinetic_energy(self):
        n3 = float(self.n) ** 6
        return 0.5 * float((self.specw[None, ...]
                            * np.abs(self.uh) ** 2).sum()) / n3

    def reynolds_stress(self):
        n3 = float(self.n) ** 6
        r = np.empty((3, 3))
        for i in range(3):
            for j in range(3):
                r[i, j] = float((self.specw
                                 * np.real(self.uh[i]
                                           * np.conj(self.uh[j]))).sum()) / n3
        return r

    def dissipation(self):
        # eps = 2 nu Omega, Omega = <omega^2>/2 = sum w k^2 |uh|^2 / (2 n^6)
        k = self.k_phys()
        k2 = (k ** 2).sum(axis=0)
        n3 = float(self.n) ** 6
        return 2.0 * self.nu * 0.5 * float(
            (self.specw[None, ...] * k2[None, ...]
             * np.abs(self.uh) ** 2).sum()) / n3

    # ------------------------------------------------------------------ #
    def _rhs(self, uh, e):
        k = self.k_phys(e)
        k2 = (k ** 2).sum(axis=0)
        k2s = k2.copy()
        k2s[0, 0, 0] = 1.0
        a = np.asarray(self.a_dir) * self.smag
        # RDT terms: -A u + 2 (k k . A u)/k^2  (A diagonal)
        au = a[:, None, None, None] * uh
        kau = (k * au).sum(axis=0) / k2s
        rhs = -au + 2.0 * k * kau[None, ...]
        if not self.linear_only:
            n = self.n
            u = np.array([np.fft.irfftn(uh[c], s=(n,) * 3, axes=(0, 1, 2)) for c in range(3)])
            nl = np.zeros_like(uh)
            for i in range(3):
                s = np.zeros((n,) * 3)
                for j in range(3):
                    dui = np.fft.irfftn(1j * k[j] * uh[i], s=(n,) * 3, axes=(0, 1, 2))
                    s += u[j] * dui
                nl[i] = np.fft.rfftn(s)
            if self.dealias:
                nl *= self.mask[None, ...]
            rhs -= nl
        rhs -= self.nu * k2[None, ...] * uh
        rhs = self.project(rhs, k)
        rhs[:, 0, 0, 0] = 0.0
        return rhs

    def step(self, dt):
        """Classical RK4 with the wavevector advanced consistently."""
        uh0, t0, e0 = self.uh, self.t, self.e
        s = self.smag

        def de(dtau):
            return s * dtau

        k1 = self._rhs(uh0, e0)
        k2 = self._rhs(uh0 + 0.5 * dt * k1, e0 + de(0.5 * dt))
        k3 = self._rhs(uh0 + 0.5 * dt * k2, e0 + de(0.5 * dt))
        k4 = self._rhs(uh0 + dt * k3, e0 + de(dt))
        self.uh = uh0 + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        self.t = t0 + dt
        self.e = e0 + de(dt)
        # keep the field exactly solenoidal wrt the new wavevector
        self.uh = self.project(self.uh, self.k_phys())
        self.uh[:, 0, 0, 0] = 0.0

    # ------------------------------------------------------------------ #
    def cfl_dt(self, cfl=0.5):
        n = self.n
        u = np.array([np.fft.irfftn(self.uh[c], s=(n,) * 3, axes=(0, 1, 2))
                      for c in range(3)])
        k = self.k_phys()
        kmax = np.array([np.abs(k[c]).max() for c in range(3)])
        umax = np.abs(u).max(axis=(1, 2, 3))
        adv = (umax * kmax).sum()
        rate = adv + (self.nu * (kmax ** 2).sum()) + abs(self.smag)
        return cfl / max(rate, 1e-12)


def vk_spectrum(amplitude=1.0, kp=4.0):
    """Passot-Pouquet-type compact spectrum peaked at kp."""
    def spec(k):
        r = k / kp
        return amplitude * r ** 4 * np.exp(-2.0 * r ** 2)
    return spec
