"""Strained-box pilot runner: decay precursor -> plane strain, with the
estimator-facing checkpoint diagnostics and the exact-RDT companion.

Protocol (cases/dns_strainedBox/SCOPING.md, decay-precursor variant):
  A. isotropic vK IC (kp, k_t = 1.5), DECAYING precursor for a fixed number
     of initial eddy turnovers (no forcing -> no forcing imprint; the
     Lee-Reynolds lineage);
  B. at switch-on measure (k_t, eps), set S = ratio * eps / k_t for the
     requested S k/eps ratio, reset e = 0, and strain to e_max, writing a
     diagnostics checkpoint (and its closed-form Cauchy-RDT companion,
     evolved from the switch-on state) at each requested e.

One realization per invocation; ensembles via the SLURM script.

Usage:
  python runner.py --n 256 --nu 2e-3 --ratio 0.8 --emax 1.0 \
      --echeck 0 0.25 0.5 0.75 1.0 --seed 1 --out /path/to/outdir
"""

from __future__ import annotations

import argparse
import copy
import os
import time

import numpy as np

from diagnostics import checkpoint
from strainbox import StrainBox, cauchy_rdt, vk_spectrum


def evolve_to(box, t_target, cfl=0.4, de_max=0.01, recompute_every=10):
    """Advance to t_target with adaptive dt (CFL + strain-increment cap)."""
    i = 0
    dt = box.cfl_dt(cfl)
    while box.t < t_target - 1e-12:
        if i % recompute_every == 0:
            dt = box.cfl_dt(cfl)
            if box.smag > 0.0:
                dt = min(dt, de_max / box.smag)
        box.step(min(dt, t_target - box.t))
        i += 1
    return i


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=256)
    ap.add_argument("--nu", type=float, default=2e-3)
    ap.add_argument("--kp", type=float, default=4.0)
    ap.add_argument("--ratio", type=float, default=0.8,
                    help="S k_t/eps at strain switch-on")
    ap.add_argument("--emax", type=float, default=1.0)
    ap.add_argument("--echeck", type=float, nargs="+",
                    default=[0.0, 0.25, 0.5, 0.75, 1.0])
    ap.add_argument("--decay-turnovers", type=float, default=1.5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=str, required=True)
    ap.add_argument("--adir", type=float, nargs=3,
                    default=[0.5, -0.5, 0.0],
                    help="strain direction cosines (traceless); "
                         "plane strain (default) or e.g. 0.25 -0.5 0.25 "
                         "for axisymmetric contraction about the line")
    ap.add_argument("--tag", type=str, default="",
                    help="extra tag for output filenames")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    tag = f"r{args.ratio:g}_s{args.seed:04d}"
    if args.tag:
        tag = args.tag + "_" + tag
    t_wall0 = time.time()

    box = StrainBox(n=args.n, nu=args.nu, smag=0.0, seed=args.seed,
                    a_dir=tuple(args.adir))
    box.init_isotropic(vk_spectrum(1.0, args.kp), kt_target=1.5)
    urms0 = np.sqrt(2.0 * box.kinetic_energy() / 3.0)
    t_pre = args.decay_turnovers / (args.kp * urms0)
    print(f"[{tag}] precursor decay to t = {t_pre:.3f}", flush=True)
    nst = evolve_to(box, t_pre)
    kt, eps = box.kinetic_energy(), box.dissipation()
    smag = args.ratio * eps / kt
    print(f"[{tag}] switch-on: kt={kt:.4f} eps={eps:.4f} "
          f"-> S={smag:.4f} ({nst} steps, "
          f"{time.time() - t_wall0:.0f} s)", flush=True)

    # reset the strain clock; labels k0 == physical k at this instant
    box.smag = smag
    box.t = 0.0
    box.e = 0.0
    uh_on = box.uh.copy()

    for e_target in sorted(args.echeck):
        if e_target > 0.0:
            nst = evolve_to(box, e_target / smag)
        path = os.path.join(args.out, f"chk_{tag}_e{e_target:g}.npz")
        d = checkpoint(box, path, extra={"seed": args.seed,
                                         "ratio": args.ratio})
        # closed-form RDT companion from the switch-on state
        rdt = copy.copy(box)
        rdt.uh = cauchy_rdt(box, e_target, uh0=uh_on)
        rdt.e, rdt.t = e_target, box.t
        checkpoint(rdt, os.path.join(args.out,
                                     f"chk_{tag}_e{e_target:g}_rdt.npz"),
                   extra={"seed": args.seed, "ratio": args.ratio,
                          "is_rdt": 1})
        print(f"[{tag}] e={e_target:g}: kt={d['kt']:.4f} "
              f"kmax_eta={d['kmax_eta']:.2f} "
              f"b={np.round(np.diag(d['R']) / np.trace(d['R']) - 1 / 3, 4)} "
              f"({time.time() - t_wall0:.0f} s)", flush=True)

    print(f"[{tag}] done in {time.time() - t_wall0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
