"""Estimator-facing diagnostics for strainbox fields (SCOPING.md section 6).

All quantities are direct spectral sums over the deforming-frame modes at
their CURRENT physical wavevectors k(t). The strain is volume-preserving
(sum a_i = 0), so the mode density in kappa-space remains exactly one mode
per unit (2 pi / L0)^3 cell and binned densities need no Jacobian.

Outputs and conventions match post/closure_bound/axisym_estimator.py:
  * phi_mn(kappa2): TWO-SIDED 1-D cross-spectral densities on the positive
    kappa2 half-axis (R_mn ~= 2 * sum phi dk2).
  * a(k2,kperp), c(k2,kperp): Batchelor-Chandrasekhar scalars from the
    azimuthally averaged spectrum tensor, via
        T = Phi_mm = 2a + c(1-x),   L = Phi_22 = (1-x)[a + c(1-x)]
        => a = T - L/(1-x),  c = (L/(1-x) - a)/(1-x),   x = k2^2/k^2.
  * Pi_rapid_exact: A_kl M_ijkl by direct quadrature (manuscript eq. 5.3).
  * m2_residue: azimuthal m=2 content of Phi_22 per (k2,kperp) bin -- the
    directly measured A1 error that section 5.4 could only bound as O(b^2).
"""

from __future__ import annotations

import numpy as np


def _mode_tensor(box):
    """Per-mode spectral tensor weights: w * Re(u_m u_n*) / n^6."""
    n6 = float(box.n) ** 6
    uh = box.uh
    w = box.specw / n6
    pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
    phi = [w * np.real(uh[m] * np.conj(uh[q])) for m, q in pairs]
    return np.array(phi)          # (6, nx, ny, nzr)


def line_spectra(box):
    """Two-sided phi_mn(kappa2) densities on the positive half-axis."""
    phi = _mode_tensor(box)
    stretch2 = np.exp(-box.a_dir[1] * box.e)
    dk2 = stretch2                                  # spacing of physical k2
    k02 = box.k0[1]
    idx = np.rint(np.abs(k02)).astype(int)
    nb = idx.max() + 1
    out = np.zeros((6, nb))
    for q in range(6):
        out[q] = np.bincount(idx.ravel(), phi[q].ravel(), minlength=nb)
    # fold: bins k>0 collect +k and -k; two-sided density = fold/2 / dk2
    out[:, 1:] *= 0.5
    kappa2 = np.arange(nb) * dk2
    return kappa2, out / dk2, dk2


def axisym_scalars(box, n_perp_bins=48):
    """Azimuthally averaged (a, c)(k2, kperp) densities plus the m=2 residue.

    Returns dict with kappa2 (n2,), kperp centers (np,), a, c (n2, np),
    m2_residue (n2, np) = |sum Phi22 e^{2 i theta}| / sum Phi22 per bin,
    and counts.
    """
    n6 = float(box.n) ** 6
    k = box.k_phys()
    k2ax = k[1]
    kperp = np.sqrt(k[0] ** 2 + k[2] ** 2)
    x = np.zeros_like(kperp)
    kk2 = k2ax ** 2 + kperp ** 2
    nz = kk2 > 0
    x[nz] = k2ax[nz] ** 2 / kk2[nz]

    w = box.specw / n6
    uh = box.uh
    trace = w * (np.abs(uh[0]) ** 2 + np.abs(uh[1]) ** 2
                 + np.abs(uh[2]) ** 2)
    lphi = w * np.abs(uh[1]) ** 2
    theta = np.arctan2(k[2], k[0])
    l_m2 = lphi * np.exp(2j * theta)

    stretch2 = np.exp(-box.a_dir[1] * box.e)
    dk2 = stretch2
    i2 = np.rint(np.abs(box.k0[1])).astype(int)
    n2 = i2.max() + 1
    kp_max = kperp.max()
    edges = np.linspace(0.0, kp_max * 1.0001, n_perp_bins + 1)
    ip = np.digitize(kperp, edges) - 1

    flat = i2 * n_perp_bins + np.clip(ip, 0, n_perp_bins - 1)
    nbins = n2 * n_perp_bins

    def acc(fld):
        return np.bincount(flat.ravel(), fld.ravel(),
                           minlength=nbins).reshape(n2, n_perp_bins)

    s_t = acc(trace)
    s_l = acc(lphi)
    s_x = acc(trace * x)
    s_m2r = acc(np.real(l_m2))
    s_m2i = acc(np.imag(l_m2))
    cnt = acc(np.ones_like(trace))

    with np.errstate(divide="ignore", invalid="ignore"):
        xb = np.where(s_t > 0, s_x / s_t, 0.0)
        # bin areas: dk2 * 2 pi kperp dkperp, x2 for the +- k2 fold
        centers = 0.5 * (edges[:-1] + edges[1:])
        area = 2.0 * dk2 * 2.0 * np.pi * centers * np.diff(edges)
        t_d = s_t / area[None, :]
        l_d = s_l / area[None, :]
        one_mx = np.clip(1.0 - xb, 1e-6, None)
        a = t_d - l_d / one_mx
        c = (l_d / one_mx - a) / one_mx
        m2 = np.where(s_l > 0,
                      np.sqrt(s_m2r ** 2 + s_m2i ** 2) / s_l, 0.0)

    return {"kappa2": np.arange(n2) * dk2, "kperp": centers,
            "a": a, "c": c, "x": xb, "m2_residue": m2, "counts": cnt}


def pi_rapid_exact(box):
    """Pi^(r)_ij = A_kl M_ijkl by direct spectral quadrature."""
    n6 = float(box.n) ** 6
    k = box.k_phys()
    k2 = (k ** 2).sum(axis=0)
    k2[0, 0, 0] = 1.0
    a_vec = np.asarray(box.a_dir) * box.smag
    uh = box.uh
    w = box.specw / n6
    # A diagonal: A_kl u_l = a_k u_k (no sum ambiguity: diag)
    au = a_vec[:, None, None, None] * uh              # (A u)_k
    pi = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            # M term 1: (k_j k_k / k^2) A_kl Phi_il -> k_j (k . A u) u_i*
            t1 = np.real((k[j] * (k * au).sum(axis=0) / k2)
                         * np.conj(uh[i]))
            t2 = np.real((k[i] * (k * au).sum(axis=0) / k2)
                         * np.conj(uh[j]))
            pi[i, j] = float((w * (t1 + t2)).sum())
    # factor 2 from the rapid Poisson source (-2 A_kl du_l/dx_k): the
    # manuscript's displayed M_ijkl (eq. 5.3) omits it -- as written it
    # yields (1/5) a1 k_t per component in the isotropic limit instead of
    # Crow's (4/5) k_t S_ij; pinned here by that exact result (see also
    # PI_PREFACTOR in post/closure_bound/axisym_estimator.py)
    return 2.0 * pi


def shell_spectrum(box):
    """E(k) on unit shells of the CURRENT physical wavenumber magnitude."""
    n6 = float(box.n) ** 6
    k = box.k_phys()
    kmag = np.sqrt((k ** 2).sum(axis=0))
    e_dens = 0.5 * (box.specw * (np.abs(box.uh) ** 2).sum(axis=0)) / n6
    shell = np.rint(kmag).astype(int)
    nb = shell.max() + 1
    ek = np.bincount(shell.ravel(), e_dens.ravel(), minlength=nb)
    return np.arange(nb), ek


def resolution_check(box):
    """kmax*eta in the worst-resolved direction (needs nu > 0)."""
    if box.nu == 0.0:
        return np.inf
    eps = box.dissipation()
    eta = (box.nu ** 3 / max(eps, 1e-300)) ** 0.25
    kmax_dir = np.array([np.abs(box.k_phys()[c]).max() for c in range(3)])
    kmax_dealiased = kmax_dir * (box.n // 3) / (box.n // 2)
    return float(kmax_dealiased.min() * eta)


def checkpoint(box, path, extra=None):
    """Write all estimator-facing diagnostics for the current state."""
    kappa2, phi, dk2 = line_spectra(box)
    ax = axisym_scalars(box)
    ks, ek = shell_spectrum(box)
    data = {"t": box.t, "e": box.e, "smag": box.smag, "nu": box.nu,
            "n": box.n, "kt": box.kinetic_energy(),
            "eps": box.dissipation(), "R": box.reynolds_stress(),
            "pi_rapid": pi_rapid_exact(box),
            "kappa2": kappa2, "phi_line": phi, "dk2": dk2,
            "shell_k": ks, "shell_E": ek,
            "kmax_eta": resolution_check(box)}
    for key, val in ax.items():
        data["ax_" + key] = val
    if extra:
        data.update(extra)
    np.savez_compressed(path, **data)
    return data
