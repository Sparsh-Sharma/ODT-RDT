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
| O8 | Minor (R1): dos Santos/Ribeiro/Piccolo refs; dissipation-estimate residual (30% discrepancy); Lee & Reynolds S-mapping justification; intro novelty overstated | OPEN | — | Introduction rewrite; dissipation estimate replaced or honestly bounded; L&R mapping discussed or the comparison demoted (R2 asks why L&R 1985 at all — consider dropping it now the DNS exists) |
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
- M0. DONE 2026-09-05 — see section 8. Lineage settled: the 25-Aug zip
  files every variant under DoNotUse/, leaving `Main_submission_V2.tex`;
  extracted to `manuscript/` (tracked).
- M1. Corrections (a)–(g) from the section5 note, folded in.
- M2. §4.2 rewrite around the three-way figure + envelope (O1, O2).
- M3. §5 rewrite to the six claims (O3, O4, O5); appendix move (R3.3).
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

Decisions (Sparsh):
- D1. Combine vs split (after M0 + a look at the length budget).
- D2. Venue.
- D3. The Kerstein-paper split (after Alan's reply).

## 7a. M0 audit findings (2026-09-05)

`manuscript/Main_submission_V2.tex` is a ~23–25 Aug post-rejection
snapshot: ~21k words (≈41 pp equivalent — unchanged from the rejected
42 pp), 15 figures, 1 table. State against the map:

- ALREADY IN V2: the strained-DNS validation section (plane strain +
  axisymmetric contraction, Reynolds-stress level — added right after the
  DNS ran); the Poisson factor-2 appears carried in the §5 source terms
  (line ~2480 — VERIFY against the note's erratum before trusting);
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
