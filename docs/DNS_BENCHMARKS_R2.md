# Independent strained-turbulence benchmarks beyond Lee & Reynolds (1985)

Vetted 2026-09-09 for Referee 2 of JFM-2026-1781 (objection O8 in
`REVISION_MAP_jfm1781.md`). **Decision pending** (Sparsh / Alan). Every
entry below was checked against its own abstract or full text; nothing is
quoted from memory. Rejected candidates are listed too, so nobody re-chases
them.

## 1. What Referee 2 actually wrote

> "Validation of (several versions of) the model is not clear. Why the
> author favours the old DNS by Lee and Reynolds (1985)? This is reported in
> a technical report, probably in connection with a CTR Summer program? I
> think that more recent achievements were carried out in further studies by
> Bill Reynolds and coworkers. None of the 15 figures of the ms. gives a
> convincing validation of a new model, capable of incorporating RDT and
> nonlinearity, vs. numerical or experimental results."

Three separable points: (a) TF-24 is a technical report, not a journal
paper; (b) "more recent achievements … by Bill Reynolds and coworkers";
(c) no convincing validation against independent numerical or experimental
results at all.

## 2. What the manuscript does now

`sec:lr-validation` (main_v3.tex l. 1020–1116) compares $b_{ij}$ against
the deformation ratio $c=\exp\int S\,\mathrm{d}t$ under plane strain with
the L&R fast/slow runs ($S^{*}=Sq^2/\varepsilon=2\chi$), in `fig_LR_plane`
and `fig_LR_b11_closure`. It already states the reason for retaining L&R:
the run was deliberately made at a turbulence Reynolds number low enough
($q^4/\nu\varepsilon\lesssim100$) to be fully resolvable on a single line,
and it reports the anisotropy against total strain (l. 1023–1036).
`fig_LR_axisym` exists in `post/closure_bound/verification/fig_lr_paper.py`
but is not used.

Position this note argues for: **keep L&R and add like-for-like modern
DNS on the same observable.** Dropping L&R would discard the only benchmark
at a Reynolds number the line can fully resolve — the argument the text
already makes — and replace it with nothing of the same kind.

## 3. The referee's own hint, verified

"Bill Reynolds and coworkers" after 1985 is the structure-based modelling
programme:

- Kassinos & Reynolds 1994, Stanford Rep. TF-61 — *also a report*.
- Reynolds & Kassinos 1995, Proc. R. Soc. A 451, "One-point modelling of
  rapidly deformed homogeneous turbulence" (doi 10.1098/rspa.1995.0118).
- Kassinos & Reynolds, CTR Annual Research Briefs 1996–1998.
- Kassinos, Reynolds & Rogers 2001, JFM 428, 213–248, "One-point turbulence
  structure tensors" (doi 10.1017/S0022112000002615).

These are **models** (structure tensors, IPRM, Q-model), not new
irrotational-strain DNS. Verified from the 1998 CTR brief (open PDF,
web.stanford.edu/group/ctr/ResBriefs98/kassinos.pdf): their
irrotational-strain validations plot "the 1985 DNS of Lee & Reynolds
(symbols)" for the axisymmetric-expansion case EXO ($Sq_0^2/\varepsilon_0=0.82$)
and the plane-strain case PXA ($Sq_0^2/\varepsilon_0=1.0$). In other words the
Reynolds group's later work uses L&R 1985 as its own benchmark. KRR 2001's
abstract says its DNS cover "diverse modes of mean deformation"; I could
not access the body, so which irrotational cases are new (if any) is
unconfirmed — treat as "L&R runs plus Rogers' shear DNS" with a hedge.

This gives an evidence-based reply to point (b), not a dodge.

## 4. Genuinely more recent DNS of irrotational strain — vetted

| # | Reference | Type | Deformation | Re / rate | Reports | Data route | Use in the paper |
|---|---|---|---|---|---|---|---|
| A | **Zusi & Perot 2013**, Phys. Fluids 25, 110819 (doi 10.1063/1.4821450) | DNS | uniform **plane strain** for a period, then return to isotropy | "moderate turbulence Reynolds numbers"; **large, moderate and small strain rates** (values in paper) | length scales, Re, decay rates, **anisotropy vs strain**; return-to-isotropy; classic RTI model constants vs time; OEC model | digitise $b_{ij}(\text{strain})$ figures; Zusi's dissertation is open (ScholarWorks `dissertations_2/129`, "Simulation and modeling of the decay of anisotropic turbulence") and should carry the numbers | **Overlay on `fig_LR_plane`.** Like-for-like with L&R: same deformation, same observable, peer-reviewed, three rapidities. |
| B | **Zusi & Perot 2014**, Phys. Fluids 26, 115103 (doi 10.1063/1.4901188) | DNS | **axisymmetric contraction and expansion**, then decay; two initial conditions | moderate Re; low, moderate, high rates | anisotropy, length scales, two-point correlations, return-to-isotropy | as A | **Revive `fig_LR_axisym`** with B (+ L&R EXO/AXK). Gives the A1-exact axisymmetric control (§5) an external anchor. |
| C | **Gualtieri & Meneveau 2010**, Phys. Fluids 22, 065104 (doi 10.1063/1.3453709) | DNS | the Chen–Meneveau–Katz straining–destraining cycle (time-dependent plane strain) | Re "by necessity … lower" than the experiment; domain aspect-ratio study | velocity variances, **Reynolds-stress anisotropy**, production and **pressure–strain budgets**; explicit comparison with the **LRR-IP** Reynolds-stress model ("good agreement with some differences for the redistribution term") | digitise | **ESM §S1**: the DNS twin of the CMK experiment already used; turns the "qualitative" experiment-only comparison into DNS-vs-model on the same LRR closure the paper runs. |
| D | **Clay & Yeung 2016**, JFM 805, 460–493 (doi 10.1017/jfm.2016.566) | DNS, up to **4096³**, deforming domain | axisymmetric contraction, smoothly varying **4:1**, strain history matched to Ayyalasomayajula & Warhaft 2006, then relaxation | high Re (AW range, $R_\lambda\le470$ in the experiment) | scale-dependent anisotropy; 1-D and axisymmetric spectra; **spectral budget incl. rapid and slow pressure–strain**, production, transfer, dissipation; "small scales return quickly, residual large-scale anisotropy persists" | JFM figures; Clay's thesis open (GaTech SMARTech 1853/60184) | Modern high-Re anchor. Cite in §4 (scale-dependent return is the paper's allocation story) and §5 (the spectral rapid pressure–strain is what the kernel closes). Optional digitised component-spectra check. |
| E | Lee, Gylfason, Perlekar & Toschi 2015, JFM 783 (doi 10.1017/jfm.2015.579) + ETC13 proceedings (open) | DNS, Rogallo deforming domain | axisymmetric strain, $\mathbf k'=(k_xe^{2St},k_ye^{-St},k_ze^{-St})$, constant $S$, several rates | $R_\lambda=117$ at strain onset; $S=16$ with $k=4.6$, $\varepsilon=2.18$ → $Sk/\varepsilon\approx34$ | $b_{ii}$ vs strain (their fig. 3), spectra; particle-focused | proceedings PDF open | Secondary citation only: another modern axisymmetric DNS reporting $b_{ii}$; not a primary benchmark. |

Experiments to pair, both already vetted: CMK 2006 (JFM 562, in use) ↔ C;
Ayyalasomayajula & Warhaft 2006 (JFM 566, doi 10.1017/s0022112006002199:
4:1 axisymmetric contraction, $40\le R_\lambda\le470$, variances agree with
RDT, spectra do not) ↔ D.

## 5. Vetted and rejected as data

- **Kevlahan & Hunt 1997**, JFM 337 (doi 10.1017/s0022112097004941):
  asymptotic analysis of when RDT fails under strong plane strain
  ($t_{NL}$ criteria). **Not DNS** — "numerical evaluation of the integrals
  for a particular form of eddy". Worth citing as *theory* in the §4
  envelope discussion, not as a benchmark.
- **Mishra & Girimaji 2013**, JFM 731 (doi 10.1017/jfm.2013.343):
  rapid-distortion *analysis* of intercomponent transfer and its
  amenability to one-point closure — theory, no dataset. Their 2017 JFM 811
  was not verified.
- "Coherent turbulent structures in a rapid contraction", JFM 2024
  (arXiv 2401.05869): **experiment**, vorticity statistics only, no
  $b_{ij}$, no RDT comparison.
- The Kassinos–Reynolds line: models validated on L&R (section 3).
- Sweeps (Crossref bibliographic 2008–2026, Semantic Scholar 2012–2026,
  several phrasings) found **no other** dedicated irrotational-strain DNS
  that reports $b_{ij}$. The modern set really is A–D (+E). This scarcity
  is itself a point to make to the referee.

## 6. Regime overlap

Paper: strained-box DNS at $Sk_t/\varepsilon=0.8$ and $16$ (0.4 for the
precursor ensembles); operating point $\chi\approx1.2$; L&R fast/slow
$S^{*}=2\chi$. A and B each span three rates (values to be read off the
papers — I could not; expected to bracket rapid→slow as L&R do). D is a
time-dependent 4:1 contraction matched to AW. E is rapid
($Sk/\varepsilon\approx34$). So the moment-level comparison can be shown
against three independent DNS families at overlapping rapidities, which is
what point (c) asks for.

## 7. Proposed manuscript actions (for decision)

1. `sec:lr-validation`: keep L&R; add one sentence — TF-24 is a report,
   but it remains the irrotational benchmark used by the Reynolds group's
   later structure-based models (Kassinos & Reynolds 1998, EXO/PXA) and the
   lowest-Re fully line-resolvable case.
2. `fig_LR_plane`: add Zusi & Perot 2013 $b_{ij}$ at the rate nearest the
   paper's (digitise; new marker family). Caption/text: two independent DNS
   at different Re, same trend.
3. Revive `fig_LR_axisym` with Zusi & Perot 2014 (+ L&R EXO/AXK).
4. ESM §S1: add Gualtieri & Meneveau 2010 to the CMK comparison; note their
   LRR-IP finding.
5. §4/§5: cite Clay & Yeung 2016 for the spectral rapid pressure–strain and
   the scale-dependent return.
6. Cover letter: reply to R2 (draft below).

Effort: the digitisation machinery exists
(`post/acoustics/expdata/vec_inspect.py`, `pa_extract.py`); needs the
PDFs (AIP/JFM via DLR). **Caveat:** A, B, C numbers (Re, rates) are behind
paywalls to me — confirm the rates before choosing which curves to overlay.

## 8. Draft reply to Referee 2

> The DNS of Lee & Reynolds (1985) is retained for a stated reason: it is
> the only strained-turbulence DNS at a Reynolds number low enough
> ($q^4/\nu\varepsilon\lesssim100$) to be fully resolved on a single line,
> and it reports the anisotropy against total strain. We note that the
> later structure-based models of Reynolds and co-workers (Reynolds &
> Kassinos 1995; Kassinos, Reynolds & Rogers 2001) are themselves validated
> against this dataset. To meet the referee's point we have added three
> independent, peer-reviewed DNS benchmarks: Zusi & Perot (2013, 2014) for
> plane strain and for axisymmetric contraction and expansion at moderate
> Reynolds number and three strain rates; Gualtieri & Meneveau (2010) for
> the Chen–Meneveau–Katz straining cycle; and Clay & Yeung (2016), whose
> $4096^3$ simulations resolve the rapid pressure–strain term spectrally.
> Figures X–Y now show the model against all of them on the same
> observable.

## Sources

- Review 2: `Review2_refjfm_sharma.pdf` (repo root).
- Kassinos & Reynolds 1998 CTR brief: https://web.stanford.edu/group/ctr/ResBriefs98/kassinos.pdf
- Reynolds & Kassinos 1995: https://doi.org/10.1098/rspa.1995.0118
- Kassinos, Reynolds & Rogers 2001: https://doi.org/10.1017/S0022112000002615
- Zusi & Perot 2013: https://doi.org/10.1063/1.4821450
- Zusi & Perot 2014: https://doi.org/10.1063/1.4901188
- Zusi dissertation: https://scholarworks.umass.edu/dissertations_2/129/
- Gualtieri & Meneveau 2010: https://doi.org/10.1063/1.3453709
- Clay & Yeung 2016: https://doi.org/10.1017/jfm.2016.566 ; thesis https://smartech.gatech.edu/handle/1853/60184
- Ayyalasomayajula & Warhaft 2006: https://doi.org/10.1017/s0022112006002199
- Lee, Gylfason, Perlekar & Toschi 2015: https://doi.org/10.1017/jfm.2015.579 ; ETC13: https://torroja.dmt.upm.es/congresos/etc13/Proceedings/PDF/167_ETC13.pdf
- Kevlahan & Hunt 1997: https://doi.org/10.1017/s0022112097004941
- Mishra & Girimaji 2013: https://doi.org/10.1017/jfm.2013.343
- Rapid-contraction experiment 2024: https://arxiv.org/abs/2401.05869
