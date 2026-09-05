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

**Scale-conditioned gate (Option A-S) DONE 2026-09-04** (commit d94b3ae:
param `anisoRejectLmax`, default 0 = bit-identical regression PASSED against
the campaign data; decks `homogeneousStrain2AS50/90`, l* = 0.05 between the
diagnostic bands; 1024 rlz each, jobs 4433518/9, ~2 min wall):

- fac 0.9 on sub-l* eddies (42% of gated candidates rejected): null except
  one marginal point (A_high at e=1: -6%, paired CI [-12%,-1%]).
- fac 0.5 on sub-l* eddies (97% rejected): A_low, A_high, transmission,
  u2^2/2kt ALL statistically unchanged.
- **Decisive extra check: high-band ENERGY is also unchanged** (paired
  median within +-25%) despite removing ~97% of small-eddy events.

**Mechanistic conclusion (the real Test-3 answer): the k2 300-800 band is
populated directly by LARGE-eddy triplet maps** — IC peak k2 ~ 50, two map
generations x3 each -> ~450 — so neither the delivery nor the relaxation of
fine-scale anisotropy runs through small eddies at this Reynolds number.
Acceptance gating of any population cannot fix transmission; the lever is
the per-map action of large eddies:
(i) Alan's **Option B** (unequal triplet-map images, middle-heavy) changes
the compression ratio per map — the direct knob;
(ii) a **scale-dependent kernel within the eddy** (kernel amplitude varying
with sub-eddy scale — the spectral analogue of the paper's g_n(x) idea; note
the allocation study already found every kernel mode kappa-rigid).

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

## 4c. Option B DONE 2026-09-05 (Alan's reply gave the spec)

Alan (2026-09-05 morning): broaden the MIDDLE image — it carries the spatial
flip (vortical idealization) and has a looser bound on attainable compression
reduction than the outers; sizes (1/6, 2/3, 1/6); keep the kernel unchanged;
outer-image boosts risk unphysical intermittency.  He also asked whether a
HiPS formulation with our results exists (it does not — clarify in reply).

Implementation (commit a7cfd6e): param `mapMidFrac` (default 1/3 =
bit-identical, regression PASSED); tripMap image volume fractions become
((1-f)/2, f, (1-f)/2); planar-only (guarded); displacement-based fillKernel()
and all tau/coefficient integrals adapt automatically.  Decks
`homogeneousStrain2B` (f=2/3) and `homogeneousStrain2B50` (f=1/2), 1024 rlz
each (jobs 4434227/4434245, ~2 min each, node-packed launcher).

**Results (`fig_optionB.png`, medians, e=3.9 unless noted):**
- **f=2/3 is the first intervention that moves anything** — at the one-point
  and energy-containing level: A_low 2.85 -> 2.50 (paired contrasts starred),
  u2^2/2kt 0.583 -> 0.559 (every e).  More kernel events per scale decade
  where the middle cascade dominates: Alan's mechanism, working as intended.
- **Fine scales unchanged; transmission NOT reduced** (0.62 -> 0.79 at
  e=3.9, because A_low improved while A_high did not): at his fractions ONE
  outer image maps the k2~50 energy peak straight into the 300-800 band
  (6 x 50 = 300) — the outer channel he called subdominant is a single-step
  feeder of unrelaxed fluid into the measurement band.
- Both bands DRAIN energy x1.4-2 (faster route to dissipation via outers).
- Intermittency (his stated worry, confirmed moderately): realizations above
  10x median high-band energy triple (13 -> 40 at e=3.9), but the largest
  single excursions are milder (max 10^2.2 vs 10^3.3) — thicker shoulder,
  tamer extremes.
- f=1/2 relaxes back to baseline except one starred A_high reduction at e=3
  (-14% paired) — the family interpolates smoothly, no sweet spot apparent.
- Structural reading: fractions sum to 1, so slowing the middle necessarily
  speeds the outers; in transmission the two nearly cancel, at one-point
  level the middle-slowing wins.
- 64-rlz "hint" (his suggestion): checked — the pseudo-success draw was
  distinguished only by absent rare high-band bursts (mean-of-ratios tail
  artifact), not a dynamical sub-ensemble worth conditioning on.

Open questions FOR Alan (in the reply): does the tighter-bound argument
imply the symmetric-outer family cannot reduce transmission at any f?
Asymmetric outers, or an eddy-type MIXTURE (mostly f=2/3 plus occasional
classic maps) as the next variant?

## 4d. Gate-A refit attempt on the 1024-rlz ensembles (2026-09-05) + mixture map

**Mixture map STAGED, not run** (commit 3156130): `mapMidFracProb` (default 1
= pure map; extra RNG draw only when a mixture is configured, so every
existing case stays bit-identical — double regression PASSED on baseline and
2B).  Per-candidate `curMidFrac` sampled in sampleEddySize, used by both
tripMaps and the LES-thirds test eddy.  Ready for Alan's answer.

**Gate-A refit on Test-3 data: NEGATIVE, and diagnostic.**  New machinery:
`dump_spectra.py` (per-k2 median + mean ensemble spectra; heavy-tail-safe)
and `gateA_fit_1024.py` (fit + Gate-B Delta-SPL, both baselines).  Verdict
(`gateA_fit_table.txt`, `fig_gateA_1024.png`): the vK-stretch family
DEGENERATES on the Test-3 ensembles — L2 -> ~0 (fit collapses to a bare
-5/3 power law), cost 45-86 in log space, Lperp/L2 in the hundreds; the
spectra are a band-limited bump + steep viscous roll-off with NO inertial
range (deck was built as a scale diagnostic, kvisc 1e-4, 8-wave IC).  Also
E1/E3 = 0.32-0.68 at e>=2: the axisymmetry assumption of the family is
badly violated by the u<->w splitting.  Do NOT use Test-3 ensembles as
Gate-A targets.

**Production case in progress:** cherry-picked tStrainOn from master
(04b66f5, partial pick of 5831da0 — allocation machinery NOT brought over);
deck `input/gateA_S1` (precursor to t=0.4, S=1, dumps e=0,0.5,1,1.5,2).
Viscosity findings: kvisc 1.44e-6 STALLS the eddy sampler (5e7 trials for
1e3 eddies, dtSmean ~4e-11, t=0.04 in 15 min); 1e-5 runs at ~2.6 h per
realization (24e6 trials by t=0.076) — feasible node-packed as ONE ~3 h
wave (raise the launcher walltime), IF the post-precursor spectrum is
fittable.  The CMK deck (caro, untracked, kvisc 1.44e-6) uses the SAME
sampler tuning (Pmax .4, Pav .02, Lp .015, Z 450; dxmin 5e-4, dxmax 1e-2) —
paper-1 runs were simply slow.

**RESOLVED same day: the stall is only the early transient** — the kvisc=1e-5
pilot completed the FULL run (precursor + e=0..2) in 14 min; the campaign
(job 4434326, walltime raised to 2:30) delivered 1024 rlz in ~35 min, 0
aborts.  `spectra_gateA_S1.npz` (median, Nu=8192) committed.

**Gate-A verdict on the production ensemble (`gateA_fit_table_gateA_S1.txt`,
`fig_gateA_gateA_S1.png`; driver now takes argv npz + E_OFF/KMAX/C0B/TAG):**
- e=0 (post-precursor): GOOD fit (cost 0.41): A0=5.3e-3, c0=-1.60, L2=11.0,
  Lperp=17.3 -> **the relaxed ODT line state is a mildly prolate vK-stretch
  (Lperp/L2 = 1.58), not isotropic** — quantified for the first time.
- e=0.5: excellent fit (cost 0.09); Gate-B: total dSPL -0.9..-3.3 dB
  (energy decay dominates), shape-only +0.9..+6.8 dB at low K_x.
- **e >= 1: the family DEGENERATES** (L2 -> 0, A0 -> 1e10, c0 slams the
  bound; widening |c0| to 8 only relocates the degeneracy).  The P2-Legendre
  h(mu) stub + rank-1 C-tensor is now DEMONSTRABLY insufficient for strained
  spectra — upgraded from "unverified stub" to data-blocked.  Consistent
  with the sec-4.2 finding (exact RDT does not rigidly translate spectra).
- **h(mu) DERIVED FROM EXACT RDT 2026-09-05** (`post/acoustics/rdt_kernel.py`,
  commit 87f3fce): plane-strain Cauchy-RDT tensor, ring-averaged about e2,
  inverted through the family's two ring equations.  h_RDT(mu;e) is NEARLY
  FLAT and positive (0.6-1.0) — the P2 stub has the WRONG SIGN for
  mu^2 < 1/3; amplitude C/A ~ e^{2e} (1.3 -> 44 by e=2); kappa-collapse in
  the inertial band; isotropy limit exact.  `axisym_family.set_h_kernel`
  hook added (default untouched, 61 tests green).
- Swapping h alone does NOT fix the axisym fits — exact RDT also violates
  the single-argument A-ansatz (fits stay degenerate with h_RDT installed).
- **Gate-A v2 WORKS** (`fit_rdt_family.py`): drop the ansatz, fit the exact
  RDT-DISTORTED vK (2 params/strain: A0, ke) to the line spectra.  Band
  [3dk, 300]: stable at ALL strains, cost 1.36-1.82, ke grows 12 -> 40 with
  e exactly as lab-frame dilatation demands.  Residual ~30% rms = the eddy
  contribution at Sk/eps ~ 0.4 (slow strain — expected).
**STALE-BINARY CORRECTION (2026-09-05 afternoon): the first gateA_S1/S20
campaigns ran a pre-tStrainOn odt.x** (the pilot script pulled but never
rebuilt) — strain + dilatation from t=0, so every number above from those
ensembles is superseded (incl. the "prolate relaxed state": that was e=0.4
data).  LESSON: any caro script that pulls a code change MUST rebuild and
re-verify before campaigning.  Clean reruns (jobs 4434552/3, precursor
verified posf0 = -0.5 at t=0.4), committed 2b11a33:

- **S1 (Sk/eps~0.4), 1024/1024**: RDT-vK fits e=0 cost 0.66 (relaxed state
  ~isotropic vK, ke 7.4; axisym fit E1/E3 = 1.002 confirms isotropy) rising
  to 1.74 at e=2 (rms ~13% -> ~21%).
- **S20 (Sk/eps~8), 626/1024** — the rapid-dilatation domainPositionToIndex
  abort takes 39% even at e<=2 (fix task queued; survivors differ at e=0:
  ke 6.2 vs 7.4 = quantified survivor bias, rapid numbers provisional).
  Fits: cost 0.58 -> 2.42 (rms ~25% at e=2).
- **KEY RESULT: the ODT-vs-exact-RDT spectral distance GROWS with Sk/eps**
  (~21% slow vs ~25% rapid at e=2) — opposite the naive rapid-limit
  expectation, exactly as the kappa-uniform rapid operator + rigid
  dilatation kinematics (sec 4.2 / Referee 3.1) predict.  First
  scale-resolved measurement of that deficiency; feeds the JFM revision
  directly.
- Axisym family on clean data: e>=1 still degenerate — structural
  conclusion stands; h(mu)-from-RDT and the RDT-vK family remain the
  working Gate-A machinery.
- NEXT: (i) fix the dilatation abort, rerun S20 full-ensemble; (ii) Gate B
  directly on exact RDT Phi_ww(kx,ky); (iii) compare the RDT-predicted
  E1/E3 split against the measured one before any non-axisymmetric
  extension.

## 5. Reply history: first reply SENT 2026-09-04 (email_alan_test3_reply.html).
   Option-B reply drafted 2026-09-05 (email_alan_optionB_reply.html).
   Original first-reply content notes below (historical):

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
- Scale-conditioned variant (our follow-up, tested): gating ONLY sub-l*
  eddies (l* = 0.05, between the bands) at 42% and even 97% rejection leaves
  fine-scale anisotropy AND fine-scale energy unchanged — the high band is
  populated directly by large-eddy maps (k2 ~ 50 -> x3 -> x3 -> ~450), so no
  acceptance-gating scheme can control transmission at this Re.  This
  sharpens his original argument: the per-map 3x compression is the whole
  story, and the counteracting relaxation must live INSIDE the large-eddy
  event.
- Therefore Option B (unequal images -> weaker compression per map) is now
  the motivated next test, or a sub-eddy scale-dependent kernel amplitude;
  both testable in the same 3-minute harness.  Ask Alan: for Option B, which
  image-size split does he want first, and does he expect the direction to
  differ between the middle-copy and outer-copy imbalance?
- The A(k) metric is model-agnostic — flow-HiPS comparison with Marten goes
  on the same plot.
- Intermittency caveat worth one sentence: rare single realizations carry
  100-2000x the median band energy, so ensemble-MEAN spectral ratios are
  non-convergent even at 1024 realizations; all quoted numbers are medians
  with bootstrap CIs (this constrains how any ODT/HiPS transmission claim
  should be measured).
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
