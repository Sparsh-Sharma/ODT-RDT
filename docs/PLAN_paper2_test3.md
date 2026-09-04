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

## 4a. CARO campaign DONE 2026-09-04 — RESULT REVERSAL

The 1024-rlz campaign ran (jobs 4431486/7, ~50 min, all COMPLETED; caro build
used the existing `~/anaconda3/envs/odt`, no micromamba needed; gate confirmed
active — runtime files report ~45% rejections).  **The 64-rlz verification in
section 3 did not survive.**

**Statistics lesson:** per-realization band energies are heavy-tailed (single
realizations carry 100–2000x the median band energy; `diag_outliers.py`).  The
ratio-of-ensemble-means A in `compare_optionA.py` does not converge even at
1024 rlz — one realization (`data_00022`, a normal completed run, i.e. genuine
ODT intermittency under strain, not a crash) alone pushed baseline A_high at
e=3.9 from ~1.75 to 0.77.  The 64-rlz numbers on both sides were tail noise.

**Robust result** (`dump_bands.py` -> `bands_*.npz` -> `robust_optionA.py`,
median-of-ratios + paired same-seed contrasts, committed as
`optionA_robust_table.txt`):

| e | baseline T = A_high/A_low | option A (fac 0.9) T | paired dlog10 A_high |
|---|---|---|---|
| 1.0 | 1.09 [1.01,1.16] | 1.02 [0.98,1.09] | -0.010 [-0.036,+0.010] |
| 2.0 | 1.00 [0.91,1.11] | 0.96 [0.90,1.03] | -0.006 [-0.035,+0.032] |
| 3.0 | 0.91 [0.80,1.06] | 0.97 [0.86,1.10] | -0.012 [-0.055,+0.043] |
| 3.9 | 0.62 [0.55,0.74] | 0.66 [0.58,0.80] | -0.005 [-0.067,+0.056] |

- **Alan's diagnosis CONFIRMED at 16x statistics**: baseline transmission ~1.0
  through e = 2–3 (imbalance reaches the small scales undiminished).
- **Option A at fac 0.9 has NO detectable effect**: every paired CI straddles
  zero (sensitivity ~+-15% on A_high — the claimed 0.96->0.54 would have been
  a -0.25 dex signal, far outside); A_low, A_high, u2^2/2kt all unchanged,
  despite the gate rejecting ~45% of accepted eddies.
- Open question: threshold too weak, or metric decoupled from A(k)?  Hence:

**fac sweep DONE 2026-09-04** (`facSweep_summary.py`, `facSweep_table.txt`,
`fig_facSweep.*`; 1024 rlz each, paired seeds; run with the node-packed
launcher `slrm_test3_node.sh` — full case in ~3 min on 8 nodes):

| fac | rejection | verdict |
|---|---|---|
| 0.9 | ~45% | indistinguishable from baseline (all bands, all e) |
| 0.8 | ~66% | indistinguishable from baseline |
| 0.7 | ~83% | indistinguishable from baseline |
| 0.5 | ~98% (~17 eddies/rlz) | A_low RISES (1.38→1.74 at e=1*), A_high rises at e≥2*, u2^2/2kt → no-eddy limit |

**Conclusion: the Option-A gate as formulated is scale-blind.**  Its
accept/reject decision barely correlates with an eddy's effect on the
spectral distribution of anisotropy; it acts as a bulk eddy-rate throttle.
At mild thresholds the surviving population reproduces the baseline
anisotropy budget exactly; at fac 0.5 the model simply drifts toward the
eddy-free (pure-LRR) limit — large-scale anisotropy grows and the fine
scales are not preferentially isotropized.  Transmission falls at fac 0.5
only because A_low rose, not because A_high dropped.  Mechanism note: as fac
drops, the accepted-candidate count RISES (279→884 per rlz) — rejected
eddies leave the field unmixed, so candidates keep firing; the gate cannot
starve itself into scale selectivity.

Natural next variants (not yet implemented): (i) **scale-conditioned gate** —
apply the rejection only to eddies smaller than a cutoff l*, directly
targeting transmission while leaving the energy-containing range alone;
(ii) Alan's **Option B** (unequal triplet-map images, middle-heavy);
(iii) a kernel-side change: make the pressure-scrambling amplitude per eddy
scale-dependent instead of gating whole eddies.

## 4b. ORIGINAL TASK (completed) — 1024-realization campaign on CARO

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

## 5. Draft reply to Alan — REWRITTEN after the 1024-rlz reversal
   (hold until the fac sweep lands; do NOT send the old 0.96->0.54 numbers)

Key content for the new draft:
- Diagnostic and campaign: A(k2)=E2/Eperp vs the model's own fixed point A=1;
  1024 paired-seed realizations per case on DLR CARO.
- His diagnosis holds with 16x statistics: median downscale transmission
  A_high/A_low = 1.09/1.00/0.91/0.62 at e = 1/2/3/3.9 — the component
  imbalance crosses ~1.7 triplet-map generations essentially undiminished
  until e ~ 4.  Per-3x-scale-step reduction ~0.8-1.0, as he estimated.
- Honest correction: our earlier 64-realization verification of Option A was
  an artifact of heavy-tailed band statistics (single realizations at 100-2000x
  the median band energy; ratio-of-means non-convergent).  With robust
  estimators and paired seeds, **fac = 0.9 produces no detectable change in
  A(k), transmission, or u2^2/2kt (sensitivity ~15%), despite rejecting ~45%
  of accepted eddies** — the surviving eddy population reproduces the same
  spectral anisotropy budget.  Interesting in itself: the gate as formulated
  (eddy-region energy-fraction metric, threshold on the post-kernel state)
  self-selects eddies whose kernels were already going to do the relaxing.
- Threshold sweep verdict (fac 0.9/0.8/0.7/0.5 = 45/66/83/98% rejection,
  1024 paired-seed rlz each): the gate is scale-blind.  0.9-0.7 are
  statistically indistinguishable from baseline in both bands; at 0.5 the
  model drifts toward the eddy-free limit (A_low up, u2^2/2kt toward 0.6006)
  without preferential fine-scale isotropization.  Interpretation: rejecting
  whole eddies modulates the rate, not the scale distribution, of
  kernel-mediated relaxation — and the unrelaxed field keeps generating
  candidates (accepted-candidate count rises 279->884 as fac drops), so the
  surviving population still executes the same cascade.
- Proposed refinement to put to Alan: condition the gate on eddy SIZE (reject
  only sub-l* eddies failing the anisotropy criterion), or make the kernel
  amplitude scale-dependent — both directly target the transmission without
  throttling the energy-containing range.  Option B (unequal images) remains
  open in the same harness.
- Option B (unequal triplet-map images) testable in the same harness; the
  A(k) metric is model-agnostic, so a flow-HiPS comparison with Marten can go
  on the same plot.
- The intermittency finding is worth a sentence on its own: rare single-eddy
  events dominate ensemble-mean fine-scale spectra under strain — relevant to
  how any anisotropy-transmission claim (ODT or HiPS) should be measured.

## 6. House rules

- Branch: `LEN_Extension` only.  Commits authored as Sparsh Sharma
  <sssparsh14@gmail.com>; no Co-Authored-By trailers per CLAUDE.md unless the
  session's harness mandates its own attribution trailer.
- Run data stays out of git (`data/**` ignored); results tables/figures and
  the scripts that regenerate them go in `post/acoustics/results_test3/`.
- Diagnostics/tests are synthetic-first: every analysis tool must prove
  itself on generated data before touching run output.
