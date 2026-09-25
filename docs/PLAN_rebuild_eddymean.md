# Rebuild: embedding the mean gradient in the eddy machinery

Branch: `eddy-mean-forcing` (off `master`; the JFM-1781 resubmission stays on
`master`). Started 2026-09-25. Driven by Alan Kerstein's objection.

## The bar Alan set

The current SC-ODT carries the mean strain as a continuous multiplicative
production `-A_ij u_j` (exact at the moment level) and applies eddy events to
the fluctuation only. Consequence: starting from laminar (`u == 0`) it stays
laminar, because the production is multiplicative and the events rearrange
zero. Alan's requirement: a correct formulation must embed the mean velocity
gradient in the eddy machinery so that, starting from `u == 0`, the mean
gradient generates turbulence — as a lab-frame ODT does, where eddy events act
on the *full* velocity including the mean profile. Acceptance test: laminar ->
turbulent transition driven by the strain, reproduced by the model.

His hypothesis: under the Rogallo transformation the eddy event itself must
change so that the mean (U_2 = -a x_2 along the line) becomes embedded in the
kernel machinery, equivalent to applying the eddy to the full velocity.

## Why this is worth doing (beyond coauthor buy-in)

The eddy-on-mean injection is **scale-selective** (energy at the eddy scales),
unlike the current wavenumber-uniform production. This is the same defect that
shows up as the fine-scale shortfall and as Referee-diagnosed "wavenumber-blind"
behaviour. If eddy-on-mean sustains the small scales on its own, the
interventions / relaxation clock of main-text sec.5-6 may become unnecessary —
a simpler, stronger paper. Risk: it is a reformulation that reopens the
validation (moments, Fig 11 spectra, LE-noise chain).

## Physics building block (STEP 1 — done, verified)

Triplet map `M` on `[y0, y0+l]`, normalized `s=(y-y0)/l`, pre-image `sigma(s)`:
`3s` on `[0,1/3]`, `2-3s` on `[1/3,2/3]` (reversed), `3s-2` on `[2/3,1]`.

Applied to the mean `U_2 = -a x_2`, the injected fluctuation
`delta = M[U_2] - U_2` is a zero-mean sawtooth with (analytic + numerically
confirmed to the digit, `post/eddymean/tripletmap_meaninject.py`):

- `<delta> = 0`            (momentum conserved; measure preserving)
- `<delta^2> = (4/27) (a l)^2`,  RMS `= 2/(3 sqrt 3) a l ~ 0.385 a l`
- energy scales exactly as `l^2`, linearly in `a`
- ~99% of the injected energy at the eddy scale (dominant at the fundamental)

So one eddy of size `l` draws `(4/27)(a l)^2` of u_2 energy from the mean, at
scale `l`. Injection is into u_2 (upwash) directly, since along the x_2 line
only U_2 has a mean gradient (U_1=U_3=0 there); kernels then redistribute.

## Remaining steps

- **STEP 2 (derivation-heavy; Opus 5).** Ensemble production rate: integrate the
  `(4/27)(a l)^2` per-eddy injection over the eddy-size density and rate
  `lambda(y0,l)` and show it recovers `P_ij = -A_ik R_kj - A_jk R_ki` in the
  existing-turbulence limit. Outcome decides whether eddy-on-mean **replaces**
  or **supplements** the current `-A_ij u_j`.
- **STEP 3 (derivation-heavy; Opus 5).** Frame transform between lab-frame ODT
  (full velocity, no compression) and the deforming frame (fluctuation +
  domain compression). Identify the extra kernel term that embeds U_2, per
  Alan's hypothesis; reconcile with the domain-straining operation (eq 2.21).
- **STEP 4.** Reimplement in the ODT solver (`code/strainbox` shares the DNS
  side; ODT changes in the fork). Run from `u == 0`, confirm generation.
- **STEP 5.** Re-validate: moments vs L&R / Zusi-Perot, spectra vs the 128^3
  strained-box DNS (Fig 11 analogue), and whether the interventions are still
  needed. Then decide the paper's new structure.

## Notes

- Keep the current forcing intact on `master`; all rebuild work here.
- Verify every derivation numerically before building on it.
- Alan is the adjudicator on ODT-orthodoxy questions.
