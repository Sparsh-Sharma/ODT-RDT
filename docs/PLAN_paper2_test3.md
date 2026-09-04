# Working context: paper 2 + Test-3 / Kerstein discussion

> Hand-off document so any Claude session (local WSL, CARO, cloud) has the full
> picture.  Updated 2026-09-04 after implementing Option A.  Everything referenced
> lives on branch `LEN_Extension`.

## 1. The two threads

**Thread A — paper 2 (leading-edge noise from strained ODT turbulence).**
Gate A: an axisymmetric anisotropy-stretched von Karman spectrum family fitted
to homogeneous-strain ODT line spectra.  Gate B: the resulting Delta-SPL of
Amiet leading-edge noise vs an isotropic baseline.  Code in `post/acoustics/`:

| module | role |
|---|---|
| `axisym_family.py` | 4-param family (A0, c0, L2, Lperp); exact isotropic vK limits |
| `fit_family.py` | params -> E_i(k2) forward map + log-space least-squares fit |
| `odt_io.py` | dmp_*.dat -> E_i(k2) spectral densities; ensemble averaging |
| `gateB_delta_spl.py` | Amiet Delta-SPL; spanwise corr. length; |L|^2 cancels |
| 61 pytest tests | all synthetic, no run data needed; `python3 -m pytest post/acoustics` |

Open Gate-A decisions (flagged as `DECISION NEEDED` in code): the h(mu)
angular kernel (currently a P2-Legendre stub) and the C-tensor form for
c0 != 0; the Delta-SPL baseline convention (`matched_isotropic` mode).
The c0 = 0 branch is exact vK and unaffected.

**Thread B — Alan Kerstein's Test-3 diagnosis (email, morning 2026-09-04).**
Test 3 = homogeneous plane strain (`input/homogeneousStrain2`: A = diag(.5,-.5,0),
S = 1, e = t, eddies ON, LRR closure, spectral band IC at 8 waves/domain,
whitened to R = (2/3)I).  Alan's argument: successive triplet maps compress
scale multiplicatively (3^N after N maps — 3 maps traverse the whole Test-3
wavenumber range) while each kernel (pressure-scrambling) event reduces
component anisotropy only weakly, so component imbalance is transported to
small scales faster than kernels can relax it.  He suggested:
- **Option A**: reject candidate eddies that do not reduce an anisotropy
  metric by a multiplicative threshold (analogous to tuning a RANS
  return-to-isotropy coefficient);
- **Option B**: unequal triplet-map images (middle larger) — direction unclear;
- flow-HiPS comparison (Marten now in the loop): HiPS swaps are even less
  scale-efficient and two identical swaps undo each other (backscatter),
  unlike triplet maps.

## 2. What we measured (this morning -> afternoon)

**Diagnostic** (`post/acoustics/scale_anisotropy.py` + `anisotropy_report.py`):
A(k2) = E2(k2)/Eperp(k2) with Eperp = (E1+E3)/2, against the **ODT-internal
isotropy reference rho_iso = 1** (`ref="equal"`): the IC gives all three
components identical line spectra and the LRR/IP kernel drives toward
component equality, so A = 1 is the model's own fixed point.  (`ref="vK"`
keeps the 3-D kinematic vK reference — rho_iso: 2 -> 3/4 — for comparison
with real turbulence; using it on ODT internals was our first mistake today.)
Decay is reported per octave (2^s) and per triplet map (3^s).

**Baseline results** (64 realizations run from source in the cloud session —
built odt.x with conda-forge cantera; ~40 s/realization;
`post/acoustics/results_test3/make_results.py`, table + figure committed):

| e | A_low (k2 30-100) | A_high (k2 300-800) | A_high/A_low | u2^2/2kt |
|---|---|---|---|---|
| 0.0 | 1.001 [0.99,1.01] | (floor) | — | 0.333 |
| 1.0 | 1.46 [1.18,1.79] | 1.09 [0.78,1.56] | 0.75 | 0.411 |
| 2.0 | 2.07 [1.71,2.49] | 1.98 [1.33,2.87] | **0.96** | 0.482 |
| 3.0 | 2.56 [1.98,3.21] | 1.76 [1.11,2.75] | 0.69 | 0.517 |
| 3.9 | 2.79 [2.30,3.33] | 1.41 [0.80,2.36] | 0.51 | 0.535 |

Bands are ~2.5-3 octaves apart (~1.7 triplet maps).  **Alan confirmed
quantitatively**: at e = 2 the imbalance reaches small scales undiminished;
per 3x scale step the anisotropy reduction factor is ~0.8-1.0.  Globally the
model is healthier: u2^2/2kt -> 0.535 vs the no-cascade LRR limit 0.6006.

## 3. Option A — implemented and verified

`src/`: `eddy::anisoReductionOK()` — after the acceptance dice roll, compute
eddy-region component energies before/after the velocity kernels analytically
on the trip-mapped copy (no domain modification):
E_after,i = E_i + c_i uRhoK_i + b_i uRhoJ_i + 0.5(c_i^2+b_i^2) rhoKK + c_i b_i rhoJK.
Metric a = |E_i/sum E - 1/3|_2; implement the eddy only if
a_after <= anisoRejectFac * a_before (regions with a_before < 0.02 pass).
Params: `LanisoReject` (default false), `anisoRejectFac` (default 0.9).
Rejection counters printed at run end.  Case `input/homogeneousStrain2A` =
baseline + gate (fac 0.9 rejects ~48% of accepted candidates).

**Verification, 64 vs 64 paired seeds** (`results_test3/compare_optionA.py`,
table/figure committed): downscale transmission A_high/A_low at
e=2: 0.96 -> 0.54; e=3: 0.69 -> 0.55; e=3.9: 0.51 -> 0.40; A_high pulled
toward 1 while A_low and u2^2/2kt are essentially unchanged — the gate
isotropizes fine scales without killing the large-scale strain response.
CIs overlap at 64 rlz, hence:

## 4. IMMEDIATE TASK — 1024-realization campaign on CARO

Sparsh wants 1024 realizations per case on CARO (DLR HPC; ssh works from the
WSL environment).  Everything is staged in `run/caro/` (see its README.md):
one-time build with micromamba conda-forge `libcantera-devel` + `yaml-cpp`
(recipe verified on Ubuntu 24.04), then

    sbatch --export=CASE=homogeneousStrain2  run/caro/slrm_test3_array.sh
    sbatch --export=CASE=homogeneousStrain2A run/caro/slrm_test3_array.sh

(128 tasks x 8 runs; set --account/--partition in the script header).  Seeds
are paired automatically (seed = 22 + shift, shift 0..1023).  Analysis:

    python3 post/acoustics/results_test3/compare_optionA.py \
        data/homogeneousStrain2 data/homogeneousStrain2A

Afterward: finalize the reply to Alan (draft below) with the 1024-run numbers,
and optionally sweep anisoRejectFac (0.95/0.8/0.7) — new input decks are a
copy of homogeneousStrain2A with one line changed.

## 5. Draft reply to Alan (update numbers after CARO, then send)

> Your scale-reduction argument prompted us to measure the competition
> directly... [full draft in the cloud-session transcript; key content:
> diagnostic definition; A_low grows to ~2.8 by e~4; transmission ~1 at e=2,
> ~0.5 by e~4 => per-map reduction ~0.8-1.0 vs the map's exact 3x compression;
> u2^2/2kt = 0.535 vs LRR 0.6006 so the deficiency is scale-local; Option A
> implemented as he proposed, threshold acts like a tunable return-to-isotropy
> coefficient — verification shows transmission at e=2 dropping 0.96 -> 0.54
> at fac=0.9 with ~48% rejections and unchanged one-point statistics;
> Option B testable in the same harness; the A(k) metric is model-agnostic,
> so a flow-HiPS comparison with Marten can go on the same plot.]

## 6. House rules

- Branch: `LEN_Extension` only.  Commits authored as Sparsh Sharma
  <sssparsh14@gmail.com>; no Co-Authored-By trailers per CLAUDE.md unless the
  session's harness mandates its own attribution trailer.
- Run data stays out of git (`data/**` ignored); results tables/figures and
  the scripts that regenerate them go in `post/acoustics/results_test3/`.
- Diagnostics/tests are synthetic-first: every analysis tool must prove
  itself on generated data before touching run output.
