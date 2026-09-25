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

## STEP 2 — ensemble production rate (done, verified; Fable 5, 2026-09-25)

**Derivation.** The triplet map is linear: `M[U+u'] = M[U] + M[u']`, so the
event on the full velocity decomposes exactly into the standard event on the
fluctuation plus the deterministic sawtooth `delta = (M-I)[U]`. The
fluctuation-energy change per event on its interval is
`int delta^2 + 2 int delta*M[u']` (the `M[u']` self-term cancels by measure
preservation). The cross term has zero ensemble mean (`<u'>=0` pointwise), so

    <DE> per event = int delta^2 = (4/27) a^2 l^3   (domain-avg: /L),
    dR22/dt|ev     = (4/27) (a^2/L) int l^3 lambda(l) dl.

Numerically confirmed (`post/eddymean/ensemble_injection.py`): per-event mean
`(4/27) a^2 l^3` at turbulence levels sigma = 0, 0.5, 2, 5 x (a l) (ratios
0.999-1.00 within SE); exact permutation map (energy preserved to machine
zero); ensemble bookkeeping over a `l^-2` size density exact to 0.06%.

**Verdict: eddy-on-mean does NOT reduce to the exact rapid production.**
`P22 = 2 a R22` is linear in `a` and proportional to `R22`; the eddy channel
is quadratic in `a` and R-independent. So it **supplements** the continuous
production; it cannot replace it. Deeper: the frame transform of the lab-frame
eddy IS `standard event + delta` (exact, from linearity of M), i.e. the
sawtooth is the eddy-level model of the SAME single exact term `u'_2 A_22`
(fluctuation-advects-mean, the only production component the line advection
can carry; the off-line components P11, P33 must stay continuous in ANY
formulation). Carrying `u'_2 A_22` continuously AND in events double-counts;
carrying it only continuously loses transition-from-laminar (current model);
carrying it only in events loses the exact rapid limit:

- ratio (eddy channel)/(exact P22) ~ O(0.1) * chi at a turbulence-set event
  rate, and steeper (~chi^2) once the mean enters the rate measure. So the
  eddy channel is negligible for chi<<1, comparable at the sustained
  equilibrium u' ~ a*l0 (there the two scalings coincide -- same physics at
  equilibrium), and would spuriously DOMINATE rapid distortion unless
  finite-time-suppressed. RDT exactness requires the eddy channel -> 0 as
  chi -> inf at fixed turbulence.

**Proposed architecture (for Alan to adjudicate).** Dose-matched split of the
one exact term: events carry the bare sawtooth `delta` (+ the mean's `A22*l`
velocity difference in the eddy rate measure, which is what makes transition
from `u'=0` possible above a Re set by `a`, `l_max`, and the viscous penalty
-- Alan's criterion); the continuous production is reduced by the fraction
`f = dR22|ev / P22` the events already deliver (all measurable on the line,
no new constant), floored at zero. Limits: `u'=0` -> pure eddy channel
(transition, Alan's bar); `chi -> inf` with ODT's existing large-eddy /
elapsed-time suppression controlling finite-time overturns -> `f -> 0`,
continuous exact production (our bar, current SC-ODT recovered). The kernel
amplitudes computed from the mapped full field automatically share the tapped
energy among components ("embedded in the kernel machinery", as Alan
guessed). General A: each component with `A_i2 != 0` gets a sawtooth with
coefficient `A_i2 * l` -- for HST (`A_12 = S`) this is exactly the standard
shear-ODT mechanism, so the formulation unifies with existing ODT practice.

## Remaining steps

- **STEP 3.** Formalise the frame transform writeup (the transform itself is
  now known: event -> event + sawtooth, mean in the rate measure, dilatation
  unchanged) and design the rapid-limit suppression + dose control precisely;
  send the note to Alan for adjudication BEFORE implementing.
- **STEP 4.** Reimplement in the ODT solver (`code/strainbox` shares the DNS
  side; ODT changes in the fork). Run from `u == 0`, confirm generation.
- **STEP 5.** Re-validate: moments vs L&R / Zusi-Perot, spectra vs the 128^3
  strained-box DNS (Fig 11 analogue), and whether the interventions are still
  needed. Then decide the paper's new structure.

## Notes

- Keep the current forcing intact on `master`; all rebuild work here.
- Verify every derivation numerically before building on it.
- Alan is the adjudicator on ODT-orthodoxy questions.
