# manuscript_compact — length-reduction workspace (JFM-2026-1781)

Parallel compression of the paper per Alan's 2026-09-16 suggestion and the
referees' length objections (R1 "shorten considerably, detail to SM";
R2 "42 pp, redundant, no clear take-home"; target ~30 pp, SOFT).

- `main_compact.tex` / `supplementary_compact.tex` — forked 2026-09-16 from
  `manuscript/main_v4_overleaf.tex` + `supplementary.tex` (identical content
  baseline, 49 pp). The live manuscript is NOT touched by work here.
- Figures are referenced from `../manuscript/Figures/` via \graphicspath
  (no duplication). Build: pdflatex+bibtex into `build/`.
- Move log:
  1. sec 6.2 intervention chronicle -> SM S2 case-by-case record;
     main text keeps the synthesis. 49 -> 48 pp (main), 6 -> 7 pp (SM).
  2. sec 2 prose trims (rate/caveat paragraph; consistency-check paragraph,
     matrix exponential inlined -- eq was unreferenced). sec 2 is otherwise
     load-bearing post-62586b2; only ~0.25 pp of true slack existed.
  3. sec 3.3 implementation verification + fig_level1a -> new SM section S5;
     main keeps compressed paragraph incl. the R1.2 LRR-definiteness statement.
  4. intro: fixed broken sentence ('six results the following aspects') and
     the red-marker CMK bullet (defects also present in the LIVE file).
  State: main 47 pp, SM 8 pp.
