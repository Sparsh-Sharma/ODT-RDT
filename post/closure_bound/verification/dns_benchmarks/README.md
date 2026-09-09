# External DNS benchmarks for the moment-level validation (R2 / O8)

Vetting and rationale: `docs/DNS_BENCHMARKS_R2.md`. Digitised data go here
as `<key>_<fig>.csv` (columns: strain-measure, value, component), one file
per curve, plus the extraction script that made them, mirroring
`../lr_extract.py` / `../lr_extracted.npz` for Lee & Reynolds.

## PDFs needed (drop them in this directory; not committed — see .gitignore)

| file name | reference | what to digitise |
|---|---|---|
| `ZusiPerot2013.pdf` | Zusi & Perot 2013, Phys. Fluids 25, 110819, doi 10.1063/1.4821450 | $b_{ij}$ (all three diagonal components) against total strain for the **large, moderate and small** strain-rate cases, plus the Reynolds number / strain-rate table |
| `GualtieriMeneveau2010.pdf` | Gualtieri & Meneveau 2010, Phys. Fluids 22, 065104, doi 10.1063/1.3453709 | velocity variances / $b_{ij}$ through the straining–destraining cycle; the LRR-IP comparison figure; the run table (Re, grid) |
| `ZusiPerot2014.pdf` | Zusi & Perot 2014, Phys. Fluids 26, 115103, doi 10.1063/1.4901188 | $b_{ij}$ against strain for axisymmetric **contraction** and **expansion** at the three rates; run table |

Clay & Yeung 2016 (JFM 805, doi 10.1017/jfm.2016.566) is cited, not
digitised, in the first pass.

## Extraction

Vector figures: `python extract_vec.py <pdf> <page>` (to be written on the
model of `post/acoustics/expdata/vec_inspect.py`); raster fallbacks use
the colour-separation route of `pa_extract.py`. Every extraction must
print a calibration check (axis tick positions round-tripped) before its
CSV is trusted, as `../fig_lr_paper.py` does for the L&R curves.

## Status (2026-09-09)

- **ZusiPerot2013 — DONE.** `extract_zp2013.py` digitises their fig. 10(a)
  (plane strain, IC3, all rates; vector, calibration residual < 1e-3).
  `zp2013_plane_high.csv` (their halved-normalisation `b_ij` relabelled to
  L&R's component convention) is overlaid as filled markers on
  `fig_LR_plane` by `../fig_lr_paper.py`.
- **ZusiPerot2014 — DONE.** `extract_zp2014.py` digitises their fig. 7
  (AXC/AXE, all rates; vector, colour→rate by the b11-extremum time).
  `zp2014_AXC_high_b11.csv` overlaid on the newly added `fig_LR_axisym`.
- **GualtieriMeneveau2010 — CITATION ONLY.** Their fig. 12 `b_ij` is a
  black-on-white RASTER of heavily overlapping hollow symbols; automated
  symbol classification was unreliable, AND the observable is the full
  straining–*destraining* cycle, which the monotonic model run does not
  reproduce (no strain reversal). So GM2010 is cited, not overlaid (ESM
  §S1). `gm2010_Sstar.csv` is the digitised fig. 1 strain history (vector,
  clean): integrating S* over the straining phase gives total deformation
  D = exp(∫S dt) ≈ 4.5, e ≈ 3.0 — recorded for reference.
