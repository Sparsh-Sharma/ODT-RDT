# Digitisation notes — far-field baseline spectra (2026-09-07)

All CSVs: two columns `f_Hz,SPL_dB`, log-spaced bin centres, median of the plotted curve
within each bin. No extrapolation beyond the plotted curves; ambiguous/overlapping bands
dropped. Each figure has a `*_check.png` re-plotting the extracted points on the original
figure in the figure's own pixel/point coordinates (visually verified).

---

## 1. Bampanis2019_AIAA.pdf — Fig 3(a), PDF page 6 (top-left of four panels)

Baseline flat plate, S_pp(dB) vs f, mic phi=0, theta=90 deg. TRUE VECTOR extraction
(`page.get_drawings()`); curves identified by stroke colour: black=(0,0,0)->19 m/s,
blue=(0,0,1)->27 m/s, red=(1,0,0)->32 m/s. Legend line samples excluded by bbox
(241,77)-(281.5,105 pt).

Calibration (page points, panel (a) axes box x 116.5-286.375, y 72.5-200.0):
- x: 185.213 pt = 10^3 Hz, 255.117 pt = 10^4 Hz (69.904 pt/decade). Cross-checked against
  minor log ticks: 136.35 pt -> 200 Hz, 148.66 pt -> 300 Hz (exact). Left spine = 104 Hz,
  right spine = 28 kHz.
- y: 200.0 pt = -40 dB, 86.667 pt = +40 dB (ticks every 20 dB); top of box = +50 dB.
  Tick values read from the 150-dpi page render.

Output (72 log bins each; bins whose median fell below -39 dB, i.e. curve clipped at the
axis floor, dropped — affects nothing except the extreme black-curve tail):
- Bampanis2019_fig3a_U19.csv: 72 pts, 116 Hz - 24.6 kHz
- Bampanis2019_fig3a_U27.csv: 72 pts, 116 Hz - 24.7 kHz
- Bampanis2019_fig3a_U32.csv: 72 pts, 116 Hz - 24.7 kHz

Uncertainty: vector centreline is essentially exact (<0.1 dB); the binned median of the
plotted narrow-band wiggles carries ~±0.6 dB (median in-bin 16-84% spread 1.2-1.3 dB) below
~8 kHz, growing to ~±2 dB in the spiky tail >15 kHz. Paper calls the data reliable only to
15 kHz — treat points above that accordingly. dB reference not stated in this paper (the
2022 JSV companion uses dB/Hz re 2e-5 Pa for the same rig).

## 2. Bampanis2022_JSV.pdf — Fig 3(b), PDF page 10 (right panel)

Baseline TIN SPL (dB/Hz re 2e-5 Pa) vs f, Theta=90 deg, Phi=0. Embedded raster
extracted at native 1333x606 px (xref 140). The three speeds are NOT colour-coded but
grey-level/weight-coded: thin black = 19 m/s, thick grey = 27 m/s, thick black = 32 m/s.
Per-column classification: grey run (110<=v<=200, not adjacent to black) -> 27 m/s;
black runs (v<90) clustered — top cluster -> 32 m/s, bottom -> 19 m/s; columns with a
single merged black cluster dropped (ambiguity rule). Vertical kc = 2pi/4pi/6pi marker
lines removed as runs taller than 60 px; legend, "(b)", kc texts and dashed arrows
excluded by rectangles.

Calibration (image px): x = 780 -> 10^2 Hz, 989 -> 10^3, 1199 -> 10^4 (209.5 px/decade;
right spine 1262 -> 2x10^4, exact). y: top frame 44 -> 60 dB, bottom 525.5 -> -20 dB
(6.02 px/dB); verified against all nine y tick-label row positions.

Output:
- Bampanis2022_fig3b_U19.csv: 56 pts, 133 Hz - 16.2 kHz
- Bampanis2022_fig3b_U27.csv: 47 pts, 110 Hz - 18.6 kHz
- Bampanis2022_fig3b_U32.csv: 57 pts, 133 Hz - 16.2 kHz

Cross-check between the two Bampanis papers (same rig, 90 deg, 32 m/s): 2019 Fig 3(a)
peak 47.5 dB @ 349 Hz vs 2022 Fig 3(b) peak 47.4 dB @ 353 Hz — consistent.

Ambiguities recorded:
- Below ~130 Hz the 19 and 32 m/s black curves merge into one cluster near the spine ->
  those columns dropped (hence 133 Hz start); the 27 m/s grey low-f stub (110-127 Hz)
  was kept; one stub-edge point at 137 Hz was dropped (median lagged the rising stub end).
- One U19 point at 10.3 kHz was dropped: it sat on the kc = 6pi vertical marker line
  inside the 19 m/s spike band (ambiguous).
- The curves are plotted with a gap ~140-200 Hz (background-subtraction dropout in the
  paper itself); no points fabricated there.
- Above ~16 kHz (19/32) the black spike bands interleave -> dropped.

Uncertainty: line half-thickness 2-3 px = ±0.4-0.5 dB; binned-median representation of the
16-Hz-bandwidth wiggles ±0.3-0.5 dB mid-band (median in-bin spread 0.5-0.9 dB), ~±1.5 dB
above ~8 kHz where the plotted curves are strongly spiky.

## 3. Narayanan2015_PoF.pdf — Fig 4, PDF page 13

U = 60 m/s flat plate. Embedded raster 767x410 px (xref 106) contains the full axes.
Legend (separate overlays on the page, read from the 150-dpi page render): solid red =
"Measured baseline flat plate" (digitised), dashed red = "Predicted baseline flat plate"
(Amiet — NOT digitised), dashed dark-orange (222,125,0) = "Background noise" (digitised).
Serrated curves (black/blue/cyan/green/magenta) excluded by colour.

Red solid vs red dashed separated by connected-component labelling: the solid curve is one
continuous component; dashes are isolated blobs. Solid then traced per-column with a
continuity tracker (max jump 8 px), which stays on the measured curve where the Amiet
dashed hump approaches it (2.5-4.5 kHz) — verified in the check overlay.

Calibration (image px): x = 100 -> 10^2 Hz, 694 -> 10^4 Hz (297 px/decade; dotted gridline
at x=397 = 10^3 confirms). y: 14 -> 75 dB, 374 -> 25 dB (7.2 px/dB, 5-dB gridlines).
Matches CASES.md stated ranges (10^2-10^4 Hz, 25-75 dB).

Output:
- Narayanan2015_fig4_baseline60.csv: 70 pts, 104 Hz - 9.1 kHz
- Narayanan2015_fig4_background60.csv: 55 pts, 107 Hz - 6.8 kHz

Notes: background curve leaves the plotted window (drops to/below the 25-dB axis floor)
beyond ~6.8 kHz -> truncated there, so the background-contaminated band is documented as
roughly f < 400 Hz (where background is within ~15 dB of the baseline) with background
known only up to 6.8 kHz. Baseline last bin at 9.1 kHz (sparse trace columns in the final
bin near 10 kHz). Observer angle for this figure is NOT STATED in the paper (CASES.md).

Uncertainty: raster 7.2 px/dB with ~2-3 px lines -> ±0.4 dB; add ±0.3 dB for the binned
median of narrow-band oscillations; overall ~±0.5-0.7 dB (background dashed similar).

### 3b. Narayanan2015 Fig 4 — Amiet prediction (red DASHED), added 2026-09-07

Same raster and calibration as section 3. The red mask splits cleanly by connected
components: solid measured = the single biggest component; the Amiet dashes = the
remaining isolated blobs (~31-44 px each). Per-column dash centres taken directly;
columns with dash runs spread > 6 px, blobs < 4 px, or dash centre within 4 px of the
solid trace dropped as ambiguous (in practice none triggered below the cutoff).

Cutoff: beyond image x ~590-600 (f ~4.5 kHz) the dashed curve merges into the thick
solid measured band (verified in a raw-pixel zoom — from ~4.5 to ~9 kHz there is only
one extra-thick red band, dashed and solid indistinguishable), and the only detached
red fragments beyond that (x >= 682, f ~9-10 kHz) sit inside the spiky multi-curve
tangle -> all dropped. Extraction therefore stops at the last clean dash (~4.4 kHz).

Output:
- Narayanan2015_fig4_amiet60.csv: 54 pts, 106 Hz - 4.2 kHz (55 log bins; med in-bin
  range 0.56 dB, max 1.11 dB — the dashed curve is smooth).

Check overlay: Narayanan2015_fig4_amiet_check.png (baseline=blue, background=lime,
Amiet=black dots); verified at zoom that black dots ride the dashes only, including
through the 2.5-4.5 kHz approach to the measured curve.

Offset (Amiet dashed minus measured solid): -5.2 dB @ 1 kHz, -5.0 dB @ 2 kHz;
at 5 kHz the dashed is merged with the measured curve in the raster (offset ~0 within
the ~±0.7 dB line width; at the last resolvable dash, 4.2-4.4 kHz, offset ~ -0.9 dB).

Uncertainty: ~±0.5 dB (±0.4 dB line centre + small binning term; no narrow-band
wiggle on this smooth prediction curve). Scripts n15_extract.py (solid/background) and
n15_amiet.py (this curve) now copied into this directory, along with the source
raster n15_p13_img0.png (embedded image xref 106 of PDF p.13; the scripts' SCR path
constant points at the original session scratchpad — repoint it here to re-run).
