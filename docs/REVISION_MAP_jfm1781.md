# Revision map — JFM-2026-1781 (rejected 18 Aug 2026)

> The working document for fixing paper 1. Objection ledger, asset inventory,
> the combine-vs-split decision, resubmission outline, and the ordered task
> list. Sources: `JFM_reject.txt` (Refs 1, 3), `Review2_refjfm_sharma.pdf`
> (Ref 2), `notes/section5_results_note.tex` (3 Sep, the quantitative core),
> `notes/allocation_results_note.tex`, and the post-3-Sep results on branch
> `LEN_Extension` (Kerstein loop + RDT-distance; see section 3).
> Manuscript sources live in `JFM_2026_ODT_RDT.zip` — note
> `Main_submission_V2.tex` there already contains post-rejection
> restructuring (strained-DNS validation sections present in its TOC);
> its exact state must be audited before writing continues (task M0).

## 0. Ground rules (from the decision letter)

- Resubmission of substantially the same paper to JFM is barred.
- ANY future JFM submission containing original material from this
  manuscript must cite JFM-2026-1781 and goes initially to the same
  associate editor (Sutanu Sarkar).
- Consequence: whatever we submit must be visibly reconstructed, and the
  cover letter must map rejection objections to changes. Venue is an open
  decision (JFM again vs JFM Rapids-to-Standard vs PoF/TCFD); all three
  referees said the §5 closure is the interesting part, so the
  reconstruction argument is strongest if §5 is the spine.

## 1. Objection ledger

Status: DONE = evidence exists and is written up in a note; DRAFT = evidence
exists, manuscript text not verified; OPEN = work or decision remains.

| # | Objection (deduplicated) | Status | Asset | Manuscript change |
|---|---|---|---|---|
| O1 | Model-vs-itself validation; no independent spectral benchmark (R1.1, R2, R3) | DONE (evidence) | 128³ strained DNS + exact-RDT companions + axisym control (`post/closure_bound/strained/`); three-way comparison (fig_threeway); **NEW: RDT-distance envelope 13→22/24% vs Sk/ε (LEN_Extension, fig_rdt_distance)** | Replace strain-on/off Figs 4–10 narrative with the three-way comparison; cite the DNS as the benchmark; add the validity-envelope figure |
| O2 | Rigid-translation "linear RDT" benchmark not derived and wrong (R3.1) | DONE — referee was right | `rdt_projection.py`: exact RDT redistributes projected upwash ×4 across the range at e=1; shape-ratio table in section5 note §5 | Correction (a): relabel eddy-free run as ODT kinematics; exact projected RDT becomes the linear reference; rewrite §4.2.1 and the abstract claim |
| O3 | §5 closure "not closed": needs a(κ₂,κ⊥), c(κ₂,κ⊥) + κ⊥ integration (R3.2) | DONE (reframed) | LOS estimator + LP bounds (los_estimator note, 21-test suite): null space stated exactly, rapid term bounded sharply; bound validity measured (e≲0.25–0.5 plane strain; exact for axisym) | §5 no longer claims a closed functional: it claims (i) exact obstruction, (ii) sharp bounds from line data, (iii) A1 measured not assumed, (iv) exact for axisymmetric distortion. **NEW: `post/acoustics/rdt_kernel.py` computes the exact (a,c) of the A1 representation from RDT — the g_n(x) kernel exercise (open item 3) is now partially in hand** |
| O4 | Most original formulation not used in the calculations (R1.2) | PART | Moment-closure reduction verified; collapsed form + A1 error evaluated on DNS spectra | State plainly which formulation generated which figure; the κ-dependent-B scoping (task C3) decides whether eq. (5.7) is exercised in ODT for the revision or stated as measured-motivation for future work |
| O5 | Rapid closure orders components wrongly under plane strain (R1.3) | DONE (diagnosed, not fixed) | Exact Π₃₃ overtaking measured (0.27 vs 0.23 at Sk/ε=16); localized in LRR b-term, not ODT; spectral c-scalar carries the ordering | New subsection: the ordering is real, sets in at rapid strain, lives in the moment closure's b-term; sensitivity statement per R1.3 |
| O6 | LE application not carried through: no Amiet mapping, no strain history (R1.4, R2) | OPEN — the combine-vs-split fork | Gate A/B machinery on LEN_Extension: RDT-vK family fits line spectra at all strains; Φ_ww(kx,ky) computable from `rdt_kernel.rdt_phi_components`; ΔSPL chain exists (gateB_delta_spl) | EITHER cut the LE claim to motivation only (split) OR add the mapping section (combine). See section 4 |
| O7 | Presentation: 42 pp, redundant, §2.1/2.2 disconnected, established-vs-new unseparated, A1 undefined, (5.8)–(5.10) underived (R2, R3.3) | OPEN | — | The restructure itself: §2 compressed to a review paragraph + SC2018 citations; A1 defined at first use; forward-map derivation to appendix; Props 1–2 demoted to consistency checks; target ≤ 30 pp |
| O8 | Minor (R1): dos Santos/Ribeiro/Piccolo refs; dissipation-estimate residual (30% discrepancy); Lee & Reynolds S-mapping justification; intro novelty overstated | OPEN | `docs/DNS_BENCHMARKS_R2.md` (candidates vetted 2026-09-09) | Introduction rewrite; dissipation estimate replaced or honestly bounded; L&R question (R2): recommendation is keep L&R + add Zusi & Perot 2013/2014, Gualtieri & Meneveau 2010, Clay & Yeung 2016 (decision pending) |
| O9 | Minor (R2): Sagaut & Cambon 2018 ch. 8 as canonical reference set; shock-turbulence intro refs off-target | OPEN | — | Adopt SC2018 as the RDT reference spine; prune intro |

## 2. What §5 (the spine) now claims — from the 3-Sep note, unchanged

1. Rapid term exact (with the Poisson factor 2 — erratum, correction b);
   single-line representation obstructed by the azimuthal null space,
   stated exactly.
2. Under A1 it collapses to the kernel form; A1 exact for axisymmetric
   distortion, measured under plane strain (residue grows with e and Sk/ε).
3. Line data bound the rapid term sharply; valid to e≈0.25–0.5 under plane
   strain — the honest statement of what a line delivers.
4. Moment closure = isotropic reduction (verified); reproduces b_ij to ~2%
   but discards the κ-dependence and the component ordering, both measured.
5. Rapid and slow terms disjoint (A=0 null + rapid-limit allocation test).
6. b₂₂(e) — the acoustically relevant quantity — linear to e=1 and
   geometry-robust.

## 3. New assets since the 3-Sep note (this week, branch LEN_Extension)

These were produced under the "paper 2 / Kerstein" label but are paper-1
material:

- **RDT-distance envelope** (`fig_rdt_distance`): rms distance of strained
  ODT line spectra from exact Cauchy-RDT-distorted vK, 13% at onset →
  ~22% (Sk/ε≈0.4) vs ~24% (Sk/ε≈8) at e=2, monotone in Sk/ε. The
  quantitative form of O1+O2's answer: ODT's rapid kinematics (κ-uniform
  operator + rigid dilatation) is where the model departs from exact
  linear theory, and the departure GROWS with strain rapidity. 1023–1024
  realizations per curve, precursor-based (tStrainOn), fixed binary.
- **Transmission diagnostic + gate verdicts** (Kerstein loop): baseline
  downscale transmission of upwash anisotropy ~1 through e≈3; acceptance
  gating (Option A, any threshold, any scale-conditioning) null; Option B
  (unequal images) moves one-point stats but not the fine scales; the high
  band is fed directly by large-eddy maps. Every measured route localizes
  the deficiency in scale-space kinematics, not relaxation efficiency.
  This is the measured backdrop for §5's g_n(x) motivation and the
  discussion of model limits; also possible material for a separate
  Kerstein co-authored paper (see section 4).
- **Exact (a,c) from RDT** (`post/acoustics/rdt_kernel.py`): ring-averaged
  exact-RDT tensor inverted through the A1 two-scalar form — h(μ) nearly
  flat and positive (the P2 guess has the wrong sign for μ²<1/3),
  amplitude ~e^{2e}. Directly feeds O3/O4 (what the κ-dependent kernel
  must look like) and kills any temptation to keep the P2 stub.
- **Dilatation abort FIXED** (f096f78 on LEN_Extension): open item 6 of the
  3-Sep note closed. High-S campaigns now run 1024/1024.
- Node-packed CARO launcher: any 1024-realization case is ~3–35 min —
  every remaining compute task below is cheap.

## 4. Combine vs split — the decision

**Option S (split, note's implicit vote):** rebuild paper 1 around §5.
Scope: strain-coupled ODT + the closure/bounds/measured-A1 story + three-way
validation + validity envelope. The LE application appears only as
motivation and one forward-looking paragraph; O6 is answered by narrowing
the claim, which R1 explicitly allows ("not sufficiently developed in
EITHER direction" → develop one fully). Paper 2 (Amiet mapping, Gate A/B,
ΔSPL) follows separately. Pros: shorter (≤30 pp reachable), strongest
reconstruction argument, fastest to submit. Cons: O6 answered by retreat;
paper 2 must then stand alone and will re-import much of paper 1.

**Option C (combine):** the above PLUS the Amiet mapping: Φ_ww(K_x,K_y)
from the exact-RDT-distorted family fitted to the line (machinery exists),
ΔSPL vs total strain and Sk/ε with the validity envelope as the error
statement. Answers all of O1–O6 head-on in one arc; clears the
"substantially different" bar most convincingly; gives the paper the
payoff R2 asked to see supported. Cons: length pressure against O7 (the
LE section costs ~5–6 pp incl. 2–3 figures); mixes a model-limits story
with an application built on the same model — must be framed as
"prediction with quantified uncertainty", which the envelope enables;
strain-history mapping to a real stagnation flow (second half of R1.4)
still deferred — must be scoped explicitly as frozen-strain ΔSPL.

**Third element — REVISED after Sparsh's 2026-09-05 call: Alan's thread
belongs to PAPER 1, not a side paper.** Alan was engaged six days after the
rejection, for the rejection; his material (Test-3 diagnosis, Option A/B
designs, allocation framework) is the "sensitivity to the known
deficiency" that R1.2/R1.3 demand. Plan of record:
- Fold the Kerstein arc into the revision COMPRESSED: one transmission
  figure + one interventions-and-verdicts table + the constructive
  endpoint (the per-map scale kinematics is the deficiency; the
  kappa-dependent kernel of eq. 5.7 is what must supply it). NOT the full
  null-result chronicle — R2's length objection stands.
- If the mixture (or a later variant) fixes transmission, the arc upgrades
  from "measured deficiency" to "deficiency -> fix" — worth waiting a few
  weeks for, not longer. Write the revision so a fix upgrades a paragraph,
  never restructures the paper (no schedule coupling to the fix-hunt).
- AUTHORSHIP: invite Alan onto the revised paper 1 — his call and
  Sparsh's; raise it in the active thread when his next reply / the
  mixture verdict lands. He may prefer acknowledgment; ask, don't assume.
- Knock-on: with Alan aboard, paper 1 is unambiguously a modeling paper —
  this tilts D1 toward Option S (LE section detaches to paper 2:
  P1 = Sharma + Kerstein model/closure/limits; P2 = Sharma Amiet
  application). A separate Kerstein methods paper (e.g. HiPS with Marten)
  becomes a CONTINGENCY for material that outgrows the revision, not the
  default.

**Coauthor recommendation (revised):** plan for Alan-in-paper-1 with the
compressed arc; hold Option C's LE section as the detachable extra it
always was (write it last; if the paper is P1+Alan, it detaches to P2).
Hard length budget ≤32 pp either way. Decision holders: Sparsh (D1, D2,
authorship invitation timing); Alan (acceptance).

## 5. Proposed outline (Option C; Option S = drop §7)

1. Introduction — distortion before the leading edge; honest novelty
   statement; dos Santos/Ribeiro/Piccolo; SC2018 as RDT spine. (3 pp)
2. Strain-coupled ODT, compact — formulation + admissibility; §2.1/2.2 of
   the old ms compressed to citations. (4 pp)
3. Verification — RDT limits, onset checks (Props 1–2 as checks). (3 pp)
4. The linear reference and the distorted spectrum — exact projected RDT
   (not rigid translation), three-way comparison ODT/RDT/DNS, the
   RDT-distance envelope vs Sk/ε. (5 pp)
5. The single-line closure (SPINE) — exact term (+factor 2), obstruction,
   A1 + measurement, LP bounds, ordering (O5), moment-closure reduction,
   disjointness; kernels/forward map in appendix. (7 pp)
6. Model limits, measured — validity envelope; κ-uniformity and the
   transmission result in one paragraph; what a κ-dependent B must supply
   (rdt_kernel result). (2 pp)
7. Leading-edge consequence [Option C only] — Amiet mapping from the
   fitted family, ΔSPL(e, Sk/ε) with envelope error bars; frozen-strain
   scoping stated. (5 pp)
8. Conclusions. (1 p)
   Appendices: forward map + kernels; numerics; IC whitening. (4 pp)

## 6. Task list (ordered)

Manuscript (master):
- fig_frozen_cost INTEGRATED 2026-09-16 (Sparsh's "money shot"; script
  `post/closure_bound/lenoise/fig_frozen_cost.py`): single panel, signed
  SPL prediction error vs total strain e — the eight fig_expval
  experiments in the laboratory window (frozen and computed coincide;
  NACA thickness cluster marked common-mode) and the exact-RDT-referenced
  band-mean errors at e=0.5/1/2 (frozen +0.6/+1.4/+3.4 dB and climbing,
  computed within 1.1–2.7 dB with measured-rms bands). Now closes Part
  III as its summary figure (fig:frozen-cost + one paragraph before the
  scope paragraph); Alan has not yet seen it — flag at next exchange.
- fig_odt_schematic REPLACED 2026-09-16 (Figure 1, sec 2.2): the
  matplotlib version is retired; the figure is now an INLINE TikZ picture
  in main_v3.tex (three panels — eddy event / space-time interleaving /
  advancement flowchart), so it inherits the JFM Times body+math fonts.
  Long descriptive sentences moved from the panels into the \caption to
  keep the panels clean and overlap-free. Canonical standalone copy of
  the TikZ (for editing outside the giant main tex) at
  `post/closure_bound/schematic/fig_odt_schematic_tikz.tex`. The
  matplotlib script and the Figures/fig_odt_schematic.pdf target were
  git-removed (superseded; recoverable from history).
- V3 SCAFFOLD DONE 2026-09-06: writing now happens in `manuscript/main_v3.tex`
  (fresh skeleton per section-5 outline, Sharma & Kerstein, provisional
  abstract limited to the six claims + envelope + hierarchical-kernel
  endpoint; every section carries TRANSPLANT markers naming its V2/notes
  sources; the two final-form M2 blocks — kinematic reference +
  exact-projected-RDT linear reference, and the strained-box DNS/three-way/
  envelope subsection — are already transplanted; compiles clean, 6 pp).
  `Main_submission_V2.tex` is FROZEN as the quarry — no further edits there.
  Remaining M-tasks execute inside main_v3.tex.
- M0. DONE 2026-09-05 — see section 8. Lineage settled: the 25-Aug zip
  files every variant under DoNotUse/, leaving `Main_submission_V2.tex`;
  extracted to `manuscript/` (tracked).
- M1. Corrections (a)–(g) from the section5 note, folded in.
- M2. RESEARCH CORE DONE 2026-09-05 (abstract/conclusions deliberately
  deferred per Sparsh): §4.2.1 relabelled to the model's kinematic
  reference + new §"linear reference for the projected spectrum" (exact
  projected RDT, fig_rdt_projection); §4.2.2 reframed (two-effects
  separation); NEW subsection sec:box-dns before §5: strained-box DNS +
  exact-RDT companions protocol, three-way table+figure (sec:threeway),
  RDT-distance envelope (sec:envelope, fig_rdt_distance). Compiles clean,
  45 pp. REMAINING in M2 scope: abstract + conclusions still carry the
  rigid-translation claim (lines ~68, ~2570-era) — do together with M4.
- M3. DONE 2026-09-06 (research core): sec:gate rebuilt in main_v3.tex — new
  intro (no closure claim), exact term + obstruction transplanted, the
  "closed computable functional" claim replaced by the collapsed-form
  status paragraph + NEW sec:bounds (LP bounds, indeterminacy table,
  T1/T2), NEW sec:a1meas (A1 measured: residue/polarisation table, O(b^2)
  superseded by measurement, bound-validity, A1 control experiment +
  fig_a1_control), NEW sec:order (Referee 1.3 measured and localised),
  lrr-reduction transplanted with measured-discards paragraph,
  doublecount + measured-disjointness paragraph, six-claims summary.
- M3b. DONE 2026-09-06: Appendix A written (app:forward, eqs. A1–A9) —
  two-scalar representation + realisability cone, forward map + T1/T2,
  kernel derivation (g_2 worked, g_1/g_3 stated; ALL THREE re-verified
  symbolically vs A_kl·M_nnkl with sympy — the §7a "VERIFY factor 2" item
  is CLOSED: eq. (5.3)-lineage conventions are consistent, no missing
  factor), LP formulation + indeterminacy + side conditions.  Lindborg1995
  added to jfm.bib.  Also fixed 3 CR-corrupted \ref{} bytes in
  sec:lrr-reduction and 2 dangling \eqref{eq:poisson} -> eq:poisson-rapid.
- M-limits. DONE 2026-09-06: sec:limits written (transmission diagnostic
  sec:transmission, interventions table tab:interventions, fig:klevels,
  endpoint sec:limits-endpoint tying scale-local relaxation to
  eq:Piclosed as the κ-dependent-operator home).  Covers the FULL 6-case
  arc incl. KS (depth tied to eddy size, best of family) and K2I3
  (iteration falsifies per-event-amplitude hypothesis); deep-strain fade
  stated as universal and left honestly undecided (artefact vs physics).
  v3 now 18 pp, compiles clean; remaining undefined refs are only the
  placeholder sections (sec:rapid-kernel, eq:rdt, sec:spec-bvc,
  sec:lr-plane).
- FIGURE STYLE (Sparsh 2026-09-06, standing rule): every plot in the
  paper must match JFM typography — same font family and PRINTED size as
  the body text.  Implemented: post/closure_bound/figstyle_jfm.py (Times/
  STIX, 9 pt base, build at true \textwidth = 32 pc = 5.31 in, never
  scale down).  Regenerated in-style: fig_klevels (new,
  post/closure_bound/limits/), fig_rdt_distance (ported from
  LEN_Extension tables), fig_pi_kernels (new analytic generator),
  fig_threeway + fig_rdt_projection (new plot-only scripts from npz),
  fig_a1_control (restyled in place).  STILL NON-COMPLIANT:
  fig_spectrum_hs__spectra/centroid (no in-repo generator — claude.ai
  era; regenerate from run data or rebuild the script) and every V2
  figure that later transplants pull in — restyle at transplant time.
- M-formulation. DONE 2026-09-06: sec:formulation transplanted into v3 —
  V2's Background folded to ONE review subsection (sec:odt-review, "a
  clear review of the JFM 2001 paper" per R2; five subsubsections
  dropped, equations kept: event/triplet/kernel/redistribution/ci +
  rate in prose); sec:rapid-slow trimmed (eq:solenoidal display -> one
  sentence); strain-forcing/rapid-kernel/domain-strain/summary at ~80%
  with Propositions 1-2 DEMOTED to consistency-check prose (R1 minor).
  §5's duplicate Poisson display removed — sec:exact-pi now recalls
  eq:poisson from §2.  New cross-link: closure-problem subsection points
  forward to sec:gate; line-realisation scope paragraph points to the
  measured wavenumber-blindness (sec:threeway).  v3 now 25 pp; dangling
  refs remaining: sec:spec-bvc, sec:lr-plane (verification + spectrum
  transplants).  Figure-style note: adopted the SoundPower jfm_rapids
  canon in figstyle_jfm (STIX 8pt, boxed axes, inward ticks, no grid,
  italic panel labels) + paper-wide colour palette — all 7 current v3
  figures regenerated.
- M-verification. DONE 2026-09-06: sec:verification transplanted (~60%)
  as "Verification, and moment-level validation" — canonical strains +
  onset merged (tab:onset kept); finite-strain with upwash necessity
  folded in (fig_rdt_components = NEW consolidated 2-panel canon figure,
  fig_upwash_amplification regenerated in canon; BOTH from
  post/closure_bound/verification/fig_rdt_moments.py, which re-verified
  every quoted number independently: onset slope to 1e-16, overtaking at
  e=3.34, u3 0.51 vs u2 0.47 at e=4, upwash peak 0.49 @ e=2.4, LRR
  0.60/IP 0.65/prod 0.98; ONE CORRECTION: axial LRR at e=4 is 0.063,
  V2 said 0.064); level1a compressed, stale caption "3%" fixed to the
  text's 0.3%; L&R compressed to one subsection (labels
  sec:lr-validation + sec:lr-plane both anchor it), R2's "why L&R"
  answered in one sentence (resolvable Re on a line; graded moment-level
  reference), S* discussion cut to one interpretive point, b11 one
  paragraph + fig_LR_b11 kept as the honesty exhibit; forward link
  finite-strain overtaking -> sec:order (now measured).  v3 31 pp; the
  ONLY dangling ref left is sec:spec-bvc (spectrum transplant).
  PENDING FIGURE RESTYLES (no generator/data locally): fig_level1a_*
  (needs a rerun of the level-1a case — also re-verify 5e-4 / 0.3%
  then), fig_LR_plane/axisym/b11_closure (need the N=1000 ensembles +
  digitised L&R points; digitised data not found in repo),
  fig_spectrum_hs__* (as before).
- FIGURE REGENERATION DONE 2026-09-07 (all 7 legacy figures now canon;
  Alan's hierarchical-reply email SENT by Sparsh 2026-09-07):
  * Deck corruption fixed: input/homogeneousStrain/input.yaml and
    run/runOneRlz.sh carried committed merge-conflict markers SINCE THE
    SUBMITTED TAG (4462f93) — the level1a deck could never have parsed.
  * Old WSL distro (all run data + N=1000 L&R ensembles + digitised
    points) is GONE (replaced by fresh Ubuntu-24.04).
  * level1a + hsA2 RERUN on caro login node (cases lvl1afig/hsA2fig,
    LEN binary, decks scp'd; dumps in scratchpad figdata).  Fresh
    verification: fractions within 7e-5 of LRR (text now "10^-4",
    tighter than V2's 5e-4), kt 0.32% ("0.3%" confirmed); kinematic
    figure: L/L0 = exp(A22 e) to 4 digits, all 3 centroids 7.03 vs
    exact 7.03 (text updated from "7.0").  New canon figures:
    fig_level1a.pdf (single 2-panel), fig_spec_kinematic.pdf (replaces
    fig_spectrum_hs__* pair, which are git-rm'd along with
    fig_level1a_trajectory*).  Scripts:
    post/closure_bound/verification/{fig_level1a,fig_spec_kinematic}.py.
  * L&R figures rebuilt from data EXTRACTED from the archived vector
    PDFs (lr_extract.py; axis-calibrated; IP/LRR curves cross-check
    against fresh recomputation to <1e-4 — calibration exact).  ODT
    curves/bands + DNS points = archived data (provenance noted in
    script); analytic curves recomputed.  fig_lr_paper.py.
  * SCIENTIFIC CORRECTION (2026-09-07, found because regeneration
    forced recomputation): V2's fig_LR_b11_closure "exact RDT" curve
    was WRONG — built from the kinematic-conservation fallacy (A_1k=0
    => R_11 conserved => b11<0).  Exact RDT (wavevector-ensemble,
    cross-checked against our own sec:finite-strain neutral-component
    result at matched deformation) gives b11 RISING to +0.114 at c=4 —
    the fast-S* DNS TRACKS exact RDT (as L&R themselves said); the
    closures (b11~0) miss the neutral-direction feed = the same
    orientation deficiency measured in sec:order.  V2's "finite-Re
    vorticity structure outside the framework" explanation is dead;
    v3 sec:lr-validation paragraph + fig:lr-b11 caption rewritten
    (with an explicit sentence retiring the wrong argument).  This
    STRENGTHENS the paper: the b11 gap is now unified with Referee
    1.3/sec:order instead of being an unexplained anomaly.
- M-spectrum. IN PROGRESS 2026-09-07: sec:spec-bvc transplanted
  (near-verbatim, reframed version); spec-Re + spec-op COMPRESSED to one
  operating-point subsubsection (eq:eps-budget kept — chi quoted from
  budget; fig_chi and fig_peak DROPPED, their content now one sentence
  each); CMK compressed (~1.5 pp: rapid limit + qualitative full-model
  + the three stated limits; hand-off to box-dns rewritten as the
  projection-matched quantitative benchmark).  ALL ensembles RERUN on
  caro (old data lost with the WSL distro): homogeneousStrainB5 (nu
  1e-5, strain on, 128 rlz), B5off (strain off baseline, 128), OP (nu
  3e-6, 16000 cells, dxmin 5e-5, 256 rlz) — decks committed in input/;
  jobs 4435510-13; postprocessor post/closure_bound/spectrum/
  ens_spectra.py (also run on caro; npz fetched).  fig_cmk_rdt DONE
  from the hsA2fig rerun: centroid 3.004 vs D=3.00; kt 1.716 (run) =
  1.719 (LRR closure, as it must) vs 1.692 exact — V2's "1.74 vs
  analytic 1.74" was sloppy, v3 states run-vs-exact honestly with the
  2% closure residual attributed.  fig_cmk_aniso now derived from the
  B5 ensemble (V2's separate N=64 ensemble had no surviving deck).
  DONE 2026-09-07 — the manuscript skeleton is COMPLETE (35 pp, ZERO
  unresolved references).  Data outcomes:
  * B5/B5off (128 rlz each, fresh): reproduce V2's claims exactly under
    V2's t=0 normalisation — full model settles 1.6-2.05, baseline
    decays to 0.64 (fig_spec_bvc canon).  fig_cmk_aniso rebuilt on a
    CLEANER construction: per-component strain-on/strain-off centroid
    ratio (cancels transient + eddy relaxation): streamwise 1.88 vs
    upwash 1.53 at e=2.2, rigid 3.0.
  * OP ensemble: caro rerun STALLED on the low-nu sampler transient
    even after the Lmin fix (e~0.13 after 3h; walltime-bound) —
    FALLBACK per Sparsh: fig_spec_op EXTRACTED from the archived vector
    PDFs (fig_spec_op_archive.py; calibration validated by the archived
    rigid line reproducing exp(e/2) to 1%; archived centroid endpoint
    2.39 = the quoted 2.4; N corrected 256->200 in text = archived
    ensemble).  chi~0.8 numbers stand on archived data.  QUEUED task:
    long-walltime OP rerun (deck input/homogeneousStrainOP with
    committed Lmin fix; also re-verify chi from eq:eps-budget then).
- SEC-6 UPDATE + APP-B + BRITISH DONE 2026-09-08: interventions table
  gains the concurrent (N=1/3) and stacked rows; fig_klevels now 8
  cases (CR3 + CRS added — Sparsh: "the concurrent thingy has a plot
  not just a table"); variants paragraph rewritten to the complete
  ledger + "resupply-limited, not implementation-limited"; Appendix B
  written (ODT configs, DNS + symmetric-C^{-1/2} whitening story, JHTDB
  nulls, statistics policy); British-English sweep (polarisation/
  realisation fixed; rest was already -ise).
- M4 (writing) LARGELY DONE 2026-09-08: introduction REWRITTEN (honest
  novelty: dosSantos2023/Ribeiro2023/Piccolo2024 + deSantana2016 cited
  as the prescribe-vs-evolve distinction per R1; SagautCambon2018 spine
  per R2; five-contribution arc paragraph; roadmap); conclusions
  WRITTEN (no rigid-translation claim; mirrors abstract); bib keys
  added and verified via publisher listings (SagautCambon2018,
  dosSantos2023 JASA 153(3):1811, Ribeiro2023 PoF 35:115112,
  Piccolo2024 PoF 36:125183); duplicate TownsendBook1976 in jfm.bib
  removed (broke bibtex).  REMAINING in M4: abstract final polish,
  LENGTH PASS (39 pp vs <=32 target — needs Sparsh's cut decisions),
  D1 (sec:le placeholder still in).
- FOUR-WAY CAPSTONE (Sparsh 2026-09-08: "I need a comparison of ODT
  standard, DNS, and my and Alan's ODT" — both models in the paper,
  head-to-head for the referees): Delta b_22(kappa_2) at e=1 on the
  threeway axes, four curves (exact RDT / DNS / strain-coupled ODT /
  ODT + adaptive-depth scale-local relaxation), BOTH Sk/eps = 0.8 (S=2
  protocol) and 16 (S=40 — this also discharges C1).  Literal
  standard ODT (no strain coupling) is a null under this protocol —
  its role is played by the strain-off baseline of sec:spec-bvc, cited
  not plotted.  Campaigns fw_S2I/S2K/S40I/S40K (fresh same-binary ISO
  baselines + KS variants, 1024 rlz, gateA_S1 precursor protocol,
  decks committed on LEN_Extension; jobs 4435645-48) RUNNING;
  analysis ready in post/closure_bound/strained/fourway.py (includes
  fresh-vs-archived ISO cross-check).  Labelling decision (variant vs
  "the model") stays with Alan's pending reply — presented in paper 1
  as the scale-local-relaxation variant.
  DONE 2026-09-08 (all 4 campaigns 1024/1024, zero bad reads):
  * Sk/eps=0.8: baseline flat (+0.07..+0.10 across bands); KS variant
    TILTS toward the DNS — fine-band +0.076 -> +0.024 vs DNS -0.046,
    roughly HALF the misallocation removed; allocation now decreases
    with wavenumber (right shape).
  * Sk/eps=16: DNS LIES ON exact RDT (rapid limit realised in the
    benchmark — new validation of the DNS+companion chain, and
    consistent with the corrected b11 story); baseline maximally wrong
    there; KS variant moves it only marginally within e<=1 — the
    resupply saturation restated on the benchmark axes.  C1 DISCHARGED.
  * CONSISTENCY FIX: archived S2_ISO (pre-dilatation-fix binary)
    disagreed with the fresh same-binary baseline in the tail-sensitive
    finest band (+0.147 -> +0.076); threeway.npz ODT rows REFRESHED
    from fw_S2I (refresh_threeway.py; archived copy kept as
    threeway_archived.npz), tab:threeway ODT row updated, TYPESw curve
    dropped from fig_threeway (different binary, unreferenced).
  * In the paper: fig:fourway + capstone paragraphs in
    sec:limits-endpoint; abstract's final sentence extended.  40 pp.
- M4. §2 compression + SC2018 + intro repair (O7, O8, O9).
- M5. [Option C] LE section from Gate machinery (O6) — written last.
- M6. Cover letter: objection→change map; cite JFM-2026-1781; venue memo.

Compute (cheap now; all on caro):
- C1. Three-way comparison at Sk/ε=16 (note open item 1; DNS data exist).
- C2. Decide on the resolved 256³ run (shared-precursor recipe exists) —
  needed only if the small-scale claims stay quantitative.
- C3. κ-dependent B scoping (note open item 3): rdt_kernel gives the
  target shape; one ODT variant + 1024-rlz campaign answers whether
  exercising eq. (5.7) is claimable in the revision or stays future work.
- C4. [Option C] ΔSPL(e, Sk/ε) from the fitted RDT-vK family — the Gate-B
  run on exact Φ_ww.

LENGTH PASS DONE 2026-09-08 (authorised): eq:triplet display -> prose;
tab:onset inlined; fig_LR_axisym dropped (numbers carried in text,
paragraph tightened); D1 DECIDED = OPTION S (sec:le placeholder deleted;
the companion-study sentence in the conclusions carries the acoustic
stage).  Result: 39 pp.  The <=32 target is NOT reachable without
structural sacrifice of load-bearing sections (candidates if forced:
sec:order or sec:doublecount compression, CMK to one figure, level1a
figure to single panel) — flagged for Sparsh; recommendation is to
argue content-per-page in the cover letter instead (39 draft-class pp
now carry DNS benchmark + bounds + measured A1 + the Kerstein arc +
four-way, vs the rejected 42 pp without any of it).

## DECISION 2026-09-08 (Sparsh, per editor guidance): ONE PAPER, THREE PARTS

The three-paper split (model / mechanism / LE application) is DEAD:
"independently these papers will lose credibility" (editor).  Everything
becomes ONE paper with three major parts + electronic supplementary
material (ESM).  Published material is trimmed to citations.

NEW ARCHITECTURE (working):
  Part I  — The model and what a line can know (existing sections 2-5,
            trimmed: standard-ODT review compressed to citations;
            CMK subsection -> ESM; level1a numerical-choices detail ->
            ESM; L&R kept condensed).
  Part II — Measured limits and the relaxation mechanism (existing
            section 6 EXTENDED with the relaxation-clock results
            (RC/RCS, currently only in the LEN plan doc) and the
            FIVE-way figure replacing the four-way; full 8-case family
            figure + intervention table detail -> ESM if length
            demands, keeping the ledger + five-way in main).
  Part III — Consequence for leading-edge noise (NEW section; the old
            Option-C content, now mandatory): the model's distorted
            upwash spectrum supplied to Amiet's response — Delta-SPL
            relative to the frozen von-Karman input, at the measured
            operating point, scoped honestly (frozen strain, no
            blocking/coherence: R1.4 second half stated as deferred,
            with the distance envelope as the uncertainty statement).
            MACHINERY EXISTS on LEN_Extension: rdt_kernel Phi_ww,
            RDT-vK family fits (fit_rdt_family), Gate-B Delta-SPL
            chain (gateB scripts) + committed clean spectra npz.
  ESM     — CMK validation; level1a numerics detail; full intervention
            family figure/table; anything else the length pass demands.

WORKING TITLE (new, to revisit): "Strain-coupled one-dimensional
turbulence for leading-edge noise: closure, measured limits, and the
scale-local relaxation mechanism".

TASK LIST (supersedes the old M-ordering):
- U1. Integrate the relaxation clock into Part II (results from
      rc/rcs tables; five-way figure into the paper; endpoint text
      updated: resupply-limited -> cured by the independent clock,
      with the kappa-operator as the remaining allocation step).
- U2. Recompute the Gate-B Delta-SPL chain on the CLEAN ensembles
      (spectra_gateA_S1/S20 npz, committed) and write Part III.
      [DONE 2026-09-09; EXTENDED 2026-09-09: absolute-SPL figure
      fig_amiet_abs added (Sparsh directive) — full Amiet chain
      (compact dipole + exact Sears, len_vsdb-audited form) at the
      sec:spec-op configuration; frozen vK baseline = e=0 member
      exactly; unit bridge validated against analytic vK (ratio
      0.79-1.14, inside fit residual) and Liepmann; distorted inflow
      lowers SPL by 5/9/14 dB at e=0.5/1/2 across the resolved band
      (<~1.8 kHz). Committed dbd722f.]
- U3. Create ESM (supplementary.tex): move CMK + level1a detail +
      full family figure; cite from main.
      [DONE 2026-09-09, commit 7754263: supplementary.tex (5 pp,
      JFM class, S-numbering) holds S1 = CMK validation (both
      figures) and S2 = full intervention family (table + 4-panel
      figure); main keeps a compact CMK summary + the prose ledger.
      level1a detail NOT moved — already compressed in the length
      pass, remainder is load-bearing.]
- U4. Trim pass: standard-ODT review to ~half via citations; abstract
      + intro + conclusions updated to the three-part claim; cover
      letter argues the one-paper reconstruction.
      [DONE 2026-09-09, commits 98f30ec + cf26b70: abstract now ends
      with the clock/five-way + the 5-20 dB Part III claim; intro =
      six contributions + explicit three-part roadmap + ESM pointer;
      conclusions paras 4-5 rewritten (fade removed by the clock,
      allocation bound stated; companion-study deferral replaced by
      the in-paper Part III result); one-paper working title
      installed in main + ESM (Sparsh to confirm); sec:odt-review
      event-rate paragraph -> citations; Cover_letter.tex rebuilt in
      the DLR template — JFM-2026-1781 cited, same-AE (Sarkar)
      request, objection-by-objection discharge, one-paper
      rationale, Kerstein as co-author.]
- U5. Length target: main text as short as honesty allows (~40 pp
      draft-class); ESM carries the rest.
      [STATUS 2026-09-09: main 41 pp + ESM 5 pp (was 42 pp rejected,
      39 pp before Part III). Close enough; only trim further if
      Sparsh wants.]
- U6. (Sparsh directive 2026-09-09) EXPERIMENTAL four-way validation:
      >=4 published far-field LE-noise cases, each with classical
      Amiet+vK, Liepmann-Amiet, standard-ODT-Amiet, fixed-ODT-Amiet
      + measured spectra overlaid; quantify which is closest.
      [IN PROGRESS. Stage 1 DONE (commit 91590d2):
       - cases: Paterson-Amiet CR-2733 (NACA0012, 40/60/90 m/s used;
         120/165 excluded, M too high for the incompressible response),
         Bampanis 2022 JSV + 2019 AIAA (ECL flat plate, 19/27/32),
         Narayanan 2015 PoF (ISVR flat plate, 60); configs + digitised
         baseline spectra in post/acoustics/expdata (CASES.md/NOTES.md,
         overlay proofs; PA fig 13 via interior-hole detection).
       - chain (post/acoustics/expval/chain.py): full non-compact
         Amiet (Bampanis 2022 eqs 2-4, L1+L2 Fresnel, mu>0.4 in all
         metric bands), 2D vK/Liepmann/ODT-representation inflows,
         Gershfeld thickness surrogate on the NACA case (all four
         chains identically), potential-flow-ellipse e_eff truncated
         at the eddy scale (PA: 0.09-0.16; flat plates: 0.03-0.05 ->
         ODT correction correctly ~vanishes there).
       - VERIFICATION: chain reproduces CR-2733's own theory curve to
         ~2-3 dB incl. the non-compactness dip; Bampanis anchor at
         27/32 m/s: bias < 0.4 dB, rms ~1.2 dB (same fidelity as
         their own Amiet validation).
       - measured so far (vk/liepmann/odt_std): ECL flat plates
         ~1 dB rms; PA overpredicted +1.8..+5.7 dB in the clean band
         (their own 1976 theory overpredicts the same data +2..+4:
         thickness+distortion); Narayanan -5 dB anomaly (observer
         angle NOT STATED in paper; their-Amiet digitisation queued
         to attribute it).
       - ENSEMBLES DONE (jobs 4437736/7): both 1024/1024, zero aborts
         (dilatation fix holds at S=20 with the clock on). KEY SIDE
         RESULT: RCS1 fit costs vs baseline - S1: 0.47-0.77 vs
         0.66-1.74; S20: 0.61-1.05 vs 0.58-2.42 (e=1: 0.97 vs 1.93;
         e=2: 0.61 vs ~2.4) - the size-capped clock moves the model
         spectra toward the exact-RDT family at every strain, both
         rapidities; npz + tables committed (5fb8683, + S20 pending
         commit).
       - DONE 2026-09-07 evening (commits 0026149, 632e850): four-way
         figures + metrics, fig_expval_paper in Part III, validation
         subsection review-hardened (test-limit framing, paired-shift
         defence, e_eff labelled estimate+band); fig_amiet_rapid added
         (exact-RDT reference, frozen error grows to +3.6 dB at e=2,
         computed inflows within 1-2.5 dB; clock gain = halved
         strain-flat uncertainty band; rotor ingestion named the
         falsifiable discriminating class); fig_dividend -> ESM S3.
         Alan emailed (validation + dividend + rotor target; sent
         2026-09-07, notes/email_alan_validation.txt); his morning
         reply: "this phase ... reached the desired outcome".]

Decisions (Sparsh):
- D1. Combine vs split (after M0 + a look at the length budget).
- D2. Venue.
- D3. The Kerstein-paper split (after Alan's reply).

## 7a. M0 audit findings (2026-09-05)

`manuscript/Main_submission_V2.tex` is a ~23–25 Aug post-rejection
snapshot: ~21k words (≈41 pp equivalent — unchanged from the rejected
42 pp), 15 figures, 1 table. State against the map:

- ALREADY IN V2: a "Validation against DNS" section that is the LEE &
  REYNOLDS 1985 comparison (b_ij level) — NOT our strained-box DNS, which
  was absent entirely (audit correction 2026-09-05); the Poisson factor-2
  appears carried in the §5 source terms
  (line ~2480 — VERIFIED 2026-09-06 via symbolic check of g_n against
  A_kl·M_nnkl: conventions consistent, see M3b);
  Sagaut & Cambon cited once.
- NOT IN V2 (everything from the 25-Aug-onward rebuild): the
  rigid-translation claim still frames the ABSTRACT (l. 68), §4.2
  (l. 1401–1492) and the CONCLUSIONS (l. 2570) — correction (a) is the
  single biggest outstanding content change, since Ref 3.1 and our own
  rdt_projection show that claim is wrong; no three-way comparison
  (correction d); no LOS-estimator/bounds reframing of §5 (O3); no
  A1-measurement / control-experiment content; no allocation verdict; no
  jackknife/realization-count reporting (f); no IC-whitening statement
  (e); no dos Santos/Ribeiro/Piccolo (O8); Lee & Reynolds still appears
  ~20 times (R2 asked why — demote now the DNS exists); length untouched
  (O7); none of this week's assets (RDT distance, transmission arc).
- Consequence for the writing order: M2 (§4.2 + abstract + conclusions
  around the correct linear reference) and M3 (§5 rebuild to the six
  claims) are the heavy lifts; M1's corrections partially reduce to them.
  The DNS section in V2 is salvageable but must be rewritten from b_ij
  level to the spectral/three-way level.

## 7b. Authorship — SETTLED

Alan Kerstein has CONFIRMED coauthorship (Sparsh, 2026-09-05). Paper 1 =
Sharma & Kerstein: the rejected manuscript's spine + the Kerstein arc
(compressed, per section 4). The paper remains open-ended by design —
ongoing research (mixture verdict, C1–C4) feeds it; structure must absorb
new results as paragraph upgrades only.

## 7. Where things live

- master: paper-1 evidence program (closure_bound, notes/*, this map).
- LEN_Extension: Kerstein loop + Gate machinery + the fixes (dilatation,
  tStrainOn partial pick). The RDT-distance and rdt_kernel assets are
  THERE — cherry-pick or merge into master when the manuscript needs the
  figures (do not regenerate by hand).
- `JFM_2026_ODT_RDT.zip` (untracked, repo root): manuscript lineage.
- CLAUDE.md now points here as the primary context; the paper-2/Kerstein
  hand-off doc remains `docs/PLAN_paper2_test3.md` on LEN_Extension.
