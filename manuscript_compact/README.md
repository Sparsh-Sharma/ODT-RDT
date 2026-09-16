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
