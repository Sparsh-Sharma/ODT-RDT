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

## Targets

- `ZusiPerot2013` → new marker family on `../fig_lr_paper.py::fig_plane`
  (`fig_LR_plane`).
- `GualtieriMeneveau2010` → ESM §S1 CMK figure (`../../spectrum/fig_cmk_*`).
- `ZusiPerot2014` → revive `../fig_lr_paper.py::fig_axisym`
  (`fig_LR_axisym`) with L&R EXO/AXK.
