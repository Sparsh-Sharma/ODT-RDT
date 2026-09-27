# HANDOFF — eddy-on-mean rebuild: the 1000-realisation re-validation campaign

Paste-in brief for the chat that holds the ODT cluster machinery (Test-3 /
Option A: the `test3n` SLURM array, the strainbox DNS, the `threeway` /
`alloc_tests` pipeline, the S2/S8 ensembles). This chat built the eddy-on-mean
rebuild but is the *paper* chat; the campaign tooling lives with you.

## 1. What the rebuild is (one paragraph)

Alan's objection: the SC-ODT model cannot generate turbulence from a laminar
(u'=0) start — the mean gradient must be embedded in the eddy machinery, not
only in a continuous production term. We implemented **eddy-on-mean**: each eddy
attempt temporarily adds the analytic mean `U_2 = A_22 (y - xc)` to `vvel`, runs
the normal eddy (rate `eddyTau`, triplet map, kernels see the full velocity),
then removes it — leaving exactly the sawtooth `(M-I)[U_2]` and putting the mean
into the eddy rate. Verified quantum: one eddy of size l injects `(4/27)(a l)^2`
into the upwash. The eddy channel is R-independent, so it **supplements** the
exact production `P_22 = 2 a R_22` (which is kept in FULL), not replaces it: at
kt->0 only the eddy fires (transition); once turbulent, production dominates.

## 2. Code state — READY TO RUN

- Branch `origin/eddymean-rebuild` (based on `LEN_Extension`, the real solver +
  clock/interventions). Tip: `3d8c802` (gate) on top of `8231647`.
- New params (all default off / inert; the model is bit-identical when off):
  - `LeddyMean` (bool, default false) — turn the mechanism on.
  - `LlaminarIC` (bool, default false) — start u=v=w=0 (transition test only).
  - `LeddyMeanKt` (double, default 1e-3) — TKE **gate** scale (see below).
- Files changed vs baseline: `src/solver.{h,cc}` (wrapper `shiftLineMean`,
  `meanGate`), `src/domain.cc` (laminar guards; production KEPT), `src/param.*`,
  `src/domaincases/domaincase_odt_homogeneousStrain.cc` (LlaminarIC).
- Build: worktree pattern on caro — `git worktree add ~/ODT-RDT-em eddymean-rebuild`,
  `conda activate odt`, `cmake -S . -B build_em -DCMAKE_BUILD_TYPE=Release ...`,
  `cmake --build build_em -j8`, `cp build_em/src/odt.x run/odt.x`. (Or just
  `bash build_odt.sh` in the worktree.)

## 3. The gate (the key refinement — and its open issue)

`meanGate()` returns `g = kref/(kref+kt)` (kt = line TKE); the injected mean is
scaled by g. Laminar (kt<<kref) -> g~1 (transition); established turbulence
(kt>>kref) -> g~0 (revert to normal ODT, no over-isotropisation). **Without the
gate the mean-in-rate over-fires eddies and over-relaxes the anisotropy** (b_22
dropped from +0.095 to +0.049 vs a DNS of +0.104 at c=1.45). With the gate,
b_22 recovers (+0.076). OPEN: `kref` is a NEW tunable constant, against this
paper's "no new constants" ethos. The clean parameter-free target is: **fire a
mean-driven eddy only when the fluctuation alone could NOT afford one** (compute
eddyTau energy with and without the mean; use fluct-only for the rate if viable,
else mean-boosted + sawtooth). Worth implementing during the campaign.

## 4. What is already established (don't re-do)

- Transition from laminar WORKS (kt: 0 -> finite, injected into R_22 first).
- Three latent bugs fixed on the branch (LRR b=R/(2kt) 0/0 at kt=0; singular
  Lyapunov at R=0; eddy-mean gate must use solver sampling `time`, not
  `mimx->time`).
- Moment b_22 (single-run, 4 seeds, L&R plane A=diag(0,-1,1), vs zp2013):
  baseline matches DNS; gated eddy ~neutral (small residual at c~1.45).
- a=2 (rapid strain) transition decays to underflow — smooth physical decay,
  low sustained level (~1e-3), NOT a numerical crash. Understand in the campaign.

## 5. THE CAMPAIGN — three decisive tests (proper statistics)

Run each with gated `LeddyMean` **on vs off**, N ~ 1024 realisations, and the
paper's own reduction (medians, centroid-normalised bands). Reuse the existing
`test3n` array + `threeway.py` + `alloc_tests.py`.

**T1 — moment anisotropy b_ij(c), confirm neutral-or-better.**
Configs: L&R plane `A=diag(0,-1,1)` and axisym; also the DNS-orientation plane
`A=diag(0.5,-0.5,0)`. Overlay Lee&Reynolds + Zusi-Perot
(`post/closure_bound/verification/dns_benchmarks/zp2013_*`). Question: gated-on
stays on top of baseline and DNS (no degradation)?

**T2 — spectral three-way (Fig 11), THE decisive test.**
Regenerate the S2/S8 ISO line-spectrum ensembles (`cases/odt_alloc/{S05,S8}_ISO`,
matched to DNS `Sk/eps = 0.8` and `16`) with gated `LeddyMean` on/off. Feed to
`post/closure_bound/strained/threeway.py` (`observables()` bins Δb_nn(κ2) in
`κ2(e)/κc(0)` bands, referenced to e=0) and compare to the 128^3 strained-box
DNS (`post/closure_bound/strained/n128/chk_r{0.8,16}_s*_e*.npz`, keys `kappa2`,
`phi_line`). QUESTION: does the scale-selective eddy injection make Δb_nn(κ2)
**wavenumber-dependent** (toward RDT/DNS) vs the baseline's wavenumber-uniform
signature? This is the Referee-1/3 objection; a clear win here justifies folding
in.

**T3 — do the interventions become unnecessary?**
With eddy-on-mean sustaining the small scales, re-run the transmission
diagnostic (`bands_homogeneousStrain2*.npz` pipeline) and check whether the
relaxation clock / sub-scale kernels (main-text sec.5-6) are still needed. If
the fine-scale imbalance A_hi stays controlled without them, the paper gets
simpler AND stronger.

## 6. PITFALLS I hit (save yourself the time)

- **Wavenumber units differ**: ODT line κ is physical (up to 1e3+); DNS κ2 is
  integer modes (1..64). MUST centroid-normalise both (the code's `centroid()`
  exists for exactly this). A raw overlay pins DNS b22 at -1/3 (garbage).
- **Spectral anisotropy is intermittent**: use MEDIANS over ~1024 realisations,
  not means; means are dominated by rare events (the paper says so).
- **Rapidity matching**: match ODT χ = Sk_t/eps to the DNS 0.8 / 16 via the
  S05/S8 configs; do not just reuse homogeneousStrain2 (χ ~ 1.2).
- **Solver vs data on different branches**: solver is on `eddymean-rebuild` (off
  LEN_Extension); the DNS data, `threeway.py`, `alloc_tests.py`, the S2/S8
  ensembles and zp csvs are on `master` (LEN_Extension is 155 behind master).
  Bring them together (cherry-pick the analysis dirs, or run the analysis from a
  master checkout pointed at the new ODT dumps).
- Runs read input from `../data/<case>/input/input.yaml` (NOT `../input/`);
  dumps to `../data/<case>/data/data_00000/dmp_*.dat`. caro python:
  `PYTHONPATH= PYTHONNOUSERSITE=1 ~/anaconda3/envs/odt/bin/python3`.

## 7. Decision criterion

"Positive" = T2 shows the gated eddy-on-mean moves the spectral anisotropy
toward the DNS (fixes the wavenumber-uniform deficiency) AND/OR T3 shows the
interventions can be dropped, with T1 confirming no moment-level degradation.
Then: draft Alan the results note (his bar + gate + these three), get his
approval, and only THEN fold into the manuscript (which stays parked on
`master`). Bring the result summary + figures back to the paper chat.
