"""Acceptance tests for the strainbox prototype.

The exact rapid-distortion reference is CLOSED-FORM (Cauchy solution for
irrotational strain): in the deforming frame, vorticity components amplify
as omega^_i(t) = exp(a_i e) omega^_i(0) at fixed mode label k0, and
u^ = i k(t) x omega^(t) / k(t)^2. This is independent of the solver's
right-hand side and time stepper, so agreement is a real validation.
"""

import numpy as np
import pytest

from strainbox import StrainBox, vk_spectrum

A_PLANE = (0.5, -0.5, 0.0)


def cauchy_rdt(box, e):
    """Exact RDT field at accumulated strain e from box's CURRENT field,
    assuming the current state is at e0 = box.e (used with e0 = 0)."""
    k0 = box.k0
    uh0 = box.uh
    om0 = 1j * np.cross(k0, uh0, axis=0)
    a = np.asarray(box.a_dir)
    om = om0 * np.exp(a * e)[:, None, None, None]
    k = box.k_phys(e)
    k2 = (k ** 2).sum(axis=0)
    k2[0, 0, 0] = 1.0
    # om = i k x u and k.u = 0  =>  k x om = -i k^2 u  =>  u = i (k x om)/k^2
    return 1j * np.cross(k, om, axis=0) / k2[None, ...]


def b_diag(box):
    r = box.reynolds_stress()
    return np.diag(r) / np.trace(r) - 1.0 / 3.0


def make_box(n=48, smag=1.0, linear=False, nu=0.0, seed=7):
    box = StrainBox(n=n, nu=nu, smag=smag, a_dir=A_PLANE,
                    linear_only=linear, seed=seed)
    box.init_isotropic(vk_spectrum(amplitude=1.0, kp=4.0), kt_target=1.5)
    return box


# --------------------------------------------------------------------------
# 0. sanity of the initial state
# --------------------------------------------------------------------------

def test_cauchy_identity_at_zero_strain():
    # u -> omega -> u round trip must be the identity at e = 0
    box = make_box(n=32)
    back = cauchy_rdt(box, e=0.0)
    err = np.sqrt(np.sum(np.abs(back - box.uh) ** 2)
                  / np.sum(np.abs(box.uh) ** 2))
    assert err < 1e-12, f"round-trip error {err:.2e}"


def test_ic_is_isotropic_and_solenoidal():
    box = make_box()
    b = b_diag(box)
    # random-phase IC with kp=4 has only ~10^2 energetic modes, so component
    # energies scatter at the few-percent level; 0.03 is the sampling floor
    assert np.all(np.abs(b) < 0.03)
    k = box.k_phys()
    div = np.abs((k * box.uh).sum(axis=0))
    scale = np.sqrt((k ** 2).sum(axis=0)) * np.abs(box.uh).max()
    assert div.max() / scale.max() < 1e-12
    assert abs(box.kinetic_energy() - 1.5) / 1.5 < 1e-12


# --------------------------------------------------------------------------
# 1. linear mode vs closed-form Cauchy RDT
# --------------------------------------------------------------------------

def test_linear_solver_matches_cauchy_exactly():
    box = make_box(n=32, smag=1.0, linear=True)
    exact = cauchy_rdt(box, e=1.0)
    nstep = 400
    dt = 1.0 / nstep
    for _ in range(nstep):
        box.step(dt)
    num, ex = box.uh, exact
    err = np.sqrt(np.sum(np.abs(num - ex) ** 2) / np.sum(np.abs(ex) ** 2))
    assert err < 1e-5, f"RDT field error {err:.2e}"


def test_rdt_energy_growth_plane_strain():
    # plane-strain RDT amplifies TKE; the Cauchy moments are the reference
    box = make_box(n=32, smag=1.0, linear=True)
    kt0 = box.kinetic_energy()
    exact = cauchy_rdt(box, e=1.0)
    box_uh_backup = box.uh.copy()
    box.uh = exact
    kt_exact = box.kinetic_energy()
    box.uh = box_uh_backup
    for _ in range(400):
        box.step(1.0 / 400)
    kt_num = box.kinetic_energy()
    assert kt_exact > kt0                       # rapid strain feeds energy
    assert abs(kt_num - kt_exact) / kt_exact < 1e-5


# --------------------------------------------------------------------------
# 2. full nonlinear solver -> RDT as S k/eps -> infinity
# --------------------------------------------------------------------------

@pytest.mark.slow
def test_nonlinear_converges_to_rdt_with_strain_rate():
    e_end = 0.5
    errs = {}
    for smag in (16.0, 128.0):
        box = make_box(n=48, smag=smag, linear=False)
        exact = cauchy_rdt(box, e=e_end)
        box_b_exact = None
        t_end = e_end / smag
        nstep = 40
        for _ in range(nstep):
            box.step(t_end / nstep)
        # compare anisotropy tensors
        bak = box.uh.copy()
        box.uh = exact
        b_ex = b_diag(box)
        box.uh = bak
        b_nu = b_diag(box)
        errs[smag] = np.abs(b_nu - b_ex).max()
    assert errs[128.0] < errs[16.0], errs
    assert errs[128.0] < 5e-3, errs


# --------------------------------------------------------------------------
# 3. conservation and structure checks on the nonlinear core
# --------------------------------------------------------------------------

def test_inviscid_unstrained_energy_conservation():
    box = make_box(n=32, smag=0.0, linear=False)
    kt0 = box.kinetic_energy()
    for _ in range(50):
        box.step(2e-3)
    drift = abs(box.kinetic_energy() - kt0) / kt0
    assert drift < 1e-8, f"energy drift {drift:.2e}"


def test_solenoidal_after_strained_nonlinear_run():
    box = make_box(n=32, smag=4.0, linear=False)
    for _ in range(30):
        box.step(2e-3)
    k = box.k_phys()
    kmag = np.sqrt((k ** 2).sum(axis=0))
    div = np.abs((k * box.uh).sum(axis=0))
    denom = (kmag * np.sqrt((np.abs(box.uh) ** 2).sum(axis=0))).max()
    assert div.max() / denom < 1e-10


def test_viscous_decay():
    box = make_box(n=32, smag=0.0, linear=False, nu=0.05)
    kt0 = box.kinetic_energy()
    for _ in range(40):
        box.step(2e-3)
    assert box.kinetic_energy() < kt0
