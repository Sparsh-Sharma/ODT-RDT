# Leading-edge (turbulence-interaction) noise experiments — baseline-case inventory

Compiled 2026-09-07 from the four PDFs in this directory. Page numbers = PDF page of the
file at hand (Narayanan/Chaitanya are author-accepted manuscripts; Bampanis files carry a
1-page HAL cover, so PDF page = printed page + 1). "NOT STATED" = not found in the paper.

---

## 1. Narayanan2015_PoF.pdf — Narayanan et al., Phys. Fluids 27, 025109 (2015)

**Facility / configuration**
- ISVR open-jet wind tunnel, anechoic chamber 8 m x 8 m x 8 m (p.8, Fig 2).
- Nozzle exit: height 0.15 m, width 0.45 m, area ratio 25:1 (p.8). Side plates at nozzle exit.
- Baseline: un-serrated FLAT PLATE, 2 mm thick, mean chord c0 = 150 mm, span 450 mm
  (p.5 states "[150 mm x 450 m]" — the "m" is a typo for mm; span L = 450 mm is confirmed p.9).
  Sharpened TE; LE tripped both sides with rough tape (p.9).
- Plate "held parallel to the jet flow direction" (p.8); AoA NOT STATED explicitly (effectively 0 deg).
- Also a NACA-65 aerofoil, same c0 and span (Sec. 3.1.2, p.13-14).

**Flow / turbulence**
- Velocities with spectra plotted: U = 20, 40, 60, 80 m/s recorded (p.9); baseline SPECTRA shown
  at U = 60 m/s (Fig 4) and U = 40 m/s (Fig 5a). Blow-down maps 30-55 m/s (Fig 10-11).
- Bi-planar grid (630 x 690 mm^2) in nozzle contraction. Hot wire at airfoil LE location,
  145 mm from nozzle exit. Tu = 2.5 %, integral length scale = 6 mm ("maintained constant"),
  obtained by fitting the measured streamwise velocity spectrum to the von Karman longitudinal
  isotropic spectrum (Pope) — p.10, Fig 3(b).
- Mean-square upwash: NOT STATED.

**Acoustics**
- 11 x 1/2" B&K 4189 mics on circular arc, R = 1.2 m from the LE, emission angles 40-140 deg
  from downstream jet axis; 10 s @ 50 kHz, window 1024 pts -> resolution 48.83 Hz (p.9).
- Quantities: SPL(f) = 10 log10(Spp(f)/pref^2), Spp = pressure PSD, pref = 20e-6 Pa (Eq 2)
  -> narrow-band PSD level in dB. PWL(f) = 10 log10(Sw(f)/Wref), Wref = 1e-12 W, sound power
  from cylindrical-radiation integration of the 40-140 deg arc, span L = 450 mm (Eqs 3-4, p.9).
- BASELINE figures:
  - Fig 4 (p.13): SPL (dB) vs Frequency (Hz), 10^2-10^4 Hz, y 25-75 dB, U = 60 m/s.
    Contains measured baseline (red solid), Amiet prediction (red dashed), background noise.
    Observer angle for this SPL spectrum: NOT STATED.
  - Fig 5(a) (p.15): PWL (dB) vs Frequency (Hz), 10^2 to 9e3 Hz, y 10-70 dB, U = 40 m/s;
    baseline flat plate AND baseline NACA-65 aerofoil (absolute levels).
- Shear-layer refraction correction: NOT STATED (none mentioned for the measurements; the
  FW-H *simulation* is corrected by (1-M^2)^2 = 0.89 for comparison, p.30).

**Vector or raster:** RASTER. Plot areas are embedded PNGs (Fig 4 plot 767x410 px; Fig 5
panels ~1150x580 px) with axis/legend text pasted as separate overlays. Readable but low-res.

**Amiet comparison:** Fig 4 (p.12-13): measured baseline flat-plate SPL vs Amiet [8] prediction:
"Good general agreement ... characteristic oscillations (interference peaks) closely matching up
to about 9 kHz"; low-frequency mismatch attributed to jet noise and grid vortex-shedding tone.

---

## 2. Chaitanya2017_JFM.pdf — Chaitanya et al., J. Fluid Mech. 818, 435-464 (2017)

**Facility / configuration**
- Same ISVR open-jet facility, anechoic chamber 8 m x 8 m x 8 m; nozzle 15 cm x 45 cm;
  aerofoil 0.15 m downstream of nozzle exit; side plates (p.7).
- Flat plates: mean chord 15 cm, span 45 cm, two 1 mm plates riveted (p.6); blunt LE, sharp TE.
- 3D aerofoil: NACA-65(12)10 (10% thickness, 1.2 lift coefficient camber), mean chord 15 cm,
  span 45 cm (p.2, p.6-7). Trip tape (SS 100, 140 um) at 16.6% chord, both sides (p.8).
- AoA = 0 deg for the acoustic results (Figs 11-13 captions); aero measurements at geometric
  -2.5 to +10 deg (effective -0.7/-1 to +2.8 deg).

**Flow / turbulence**
- U = 20, 40, 60, 80 m/s (p.8). Most spectra shown at 60 m/s; Fig 10 baseline PWL at 20 and 60 m/s.
- Grid 1: Tu = 2.5 %, streamwise integral length scale 7.5 mm; Grid 2 (holes 40 mm): Tu = 4.5 %,
  streamwise scale ~13.5 mm. vK fit to measured streamwise spectrum at 145 mm downstream of
  nozzle exit; TRANSVERSE scale Lambda_t = streamwise/2 = 3.75 mm and 6.75 mm (p.8-9, Fig 3).
- Mean-square upwash: NOT STATED.

**Acoustics**
- 11 x 1/2" B&K 4189 mics, R = 1.2 m from mid-span LE, 40-140 deg; 10 s @ 50 kHz,
  resolution 48.83 Hz; SPL/PWL computed as in Narayanan et al. 2015 (p.8).
- BASELINE absolute spectra: ONLY Figure 10 (p.20): PWL (dB) vs non-dimensional frequency
  fc0/U (x 10^-1 to ~2x10^2, y -10 to 70 dB), baseline NACA65 total noise + self-noise at
  U = 20 and 60 m/s. NO absolute baseline flat-plate spectrum vs Hz appears in this paper —
  all other acoustic figures are reduction spectra (DeltaPWL, DeltaOAPWL) or per-valley power.
- Shear-layer refraction correction: NOT STATED.

**Vector or raster:** RASTER (one embedded PNG per figure, e.g. Fig 10 = 819x569 px).

**Amiet comparison:** none of their own measurements vs Amiet in this paper. Text (p.5) notes
Lyu et al. (2016) serration model (iterative Amiet approach) gives "predictions ... in close
agreement with the experimental data". No baseline-vs-Amiet figure.

---

## 3. Bampanis2019_AIAA.pdf — Bampanis, Roger, Ragni, Avallone & Teruna, AIAA 2019-2741

**Facility / configuration**
- ECL (Ecole Centrale de Lyon) low-speed anechoic open-jet tunnel, room 4 m x 5 m x 6 m;
  rectangular nozzle, vertical outlet cross-section 15 cm x 30 cm (PDF p.3).
- FLAT-PLATE airfoils: NACA symmetric 4-digit profile split at max thickness with flat mid part;
  thickness 3 mm, mean chord c0 = 10 cm, wetted span 30 cm (= nozzle height). Baseline =
  straight LE; two serrated variants (h/c0 = 0.07 approx./0.1; deepest h/c0 = 0.167,
  inclination 76 deg) (PDF p.3-4). Held on narrow supports (not end plates).
- AoA = 0 deg, fixed ("the zero angle of attack cannot be varied with the present setup", PDF p.4).

**Flow / turbulence**
- Baseline far-field spectra plotted at U = 19, 27 and 32 m/s (Fig 3a).
- Grid upstream of contraction (bars 8 mm x 2 mm, mesh 5 cm): Tu = 4.5 %, integral length
  scale = 9 mm, from fitting a model von Karman spectrum to the measured streamwise velocity
  spectrum (single hot wire at the LE location, airfoil absent) (PDF p.3).
- Mean-square upwash: not measured directly; model input sqrt(u2bar) = TI*U/100 (Eq 2, PDF p.9).

**Acoustics**
- Rotating vertical arc of 6 x 1/2" B&K 4189 mics; R = 1.25 m from LE center; elevation
  phi = 0-75 deg (steps 15 deg), meridian angle theta = 20-110 deg from streamwise direction;
  51.2 kHz, 30 x 1 s averages, resolution 1 Hz (PDF p.3).
- BASELINE spectra: Fig 3(a) (PDF p.6): "S_pp (dB)" vs Frequency (Hz), log axis roughly
  2x10^2 to 2.5x10^4 Hz (data called reliable to 15 kHz), y -40 to +50 dB; U = 19/27/32 m/s;
  microphone phi = 0 deg (mid-span), theta = 90 deg. PSD of acoustic pressure with background
  noise subtracted; plotted with bandwidth reduced to 16 Hz. dB reference NOT STATED on the
  axis (the companion JSV 2022 paper uses dB/Hz re 2e-5 Pa for the same setup).
  Also Fig 5(b) (PDF p.7): total noise vs TE-noise spectra, 32 m/s, 90 deg, mid-span.
- Shear-layer refraction: "corrections were not implemented for the refraction of sound ...
  These corrections are negligible for the present experiment" (PDF p.10).

**Vector or raster:** VECTOR — all plot pages are pure vector drawings (text/curves stay sharp).

**Amiet comparison:** Fig 6 (PDF p.8): measured 3D directivity (0.5-7.2 kHz) vs Amiet TIN model
(blue lines), baseline only: "good performance in all directions", but "Amiet's model
underestimates the sound radiation at low frequencies" (PDF p.10).

---

## 4. Bampanis2022_JSV.pdf — Bampanis, Roger & Moreau, J. Sound Vib. 520, 116635 (2022)

**Facility / configuration**
- ECL anechoic open-jet: 5 m (flow) x 6 m (lateral) x 4 m (height); cut-off < 100 Hz.
- Nozzle: contraction 2:1, from 30 x 30 cm to vertical rectangular 30 cm x 15 cm (PDF p.7).
- Baseline: FLAT PLATE 3 mm thick (e/c = 3%), chord c = 10 cm, span L = 30 cm (= nozzle
  height); rounded LE + sharp TE from a NACA-0004 profile cut at max thickness (PDF p.7, Fig 1b).
  (NOT a NACA-0012 — that airfoil belongs to their porous-airfoil companion studies.)
- Narrow supports (no end plates); AoA = 0 deg; Re_c = 1.3e5-2.1e5 (U0 = 19-32 m/s).

**Flow / turbulence**
- Baseline TIN spectra plotted at U0 = 19, 27, 32 m/s (Fig 3b); all serration/directivity/Amiet
  results shown at 32 m/s only ("similar conclusions ... for the other speeds", PDF p.11).
- Same grid (mesh 5 cm, flat bars 8 mm). Hot wire 10 cm downstream of nozzle exit at the LE
  location, airfoil removed: Lambda = 9 mm and u_rms/U0 = 4.5 %, determined by tuning the von
  Karman model spectrum on the measured streamwise-velocity PSD, with an added high-frequency
  viscous correction factor exp[-0.0008 (k1/ke)^2] (PDF p.8, Fig 2).
- Mean-square upwash: via isotropic vK Phi_ww with u_rms above (PDF p.8); not separately measured.

**Acoustics**
- Arc antenna of 6 x 1/2" B&K 4189 mics, Rm = 1.25 m; azimuth Theta accessible -110 to +110 deg
  (results at 30/60/90/110 deg), elevation Phi = 0-75 deg by 15 deg; 51.2 kHz, 30 x 1 s,
  1 Hz resolution reduced to 16 Hz bandwidth after background subtraction (PDF p.7, p.9).
- Quantity: SPL (dB/Hz, ref. 2 10^-5 Pa) — sound-pressure PSD level. TIN isolated by subtracting
  background noise, then the separately measured self-noise (grid removed, tripped) (PDF p.9).
- BASELINE figures:
  - Fig 3 (PDF p.10): SPL (dB/Hz re 2e-5 Pa) vs frequency (Hz), 10^2 to ~2x10^4 Hz, y -20 to
    60 dB; (a) subtraction demo at 32 m/s; (b) baseline TIN at 19/27/32 m/s; Theta = 90 deg,
    Phi = 0 deg.
  - Figs 6-7 (PDF p.13-14): baseline vs serrated TIN spectra at Theta = 30/60/90/110 deg,
    Phi = 0/45/60 deg, 32 m/s (pairs vertically shifted for clarity).
  - Fig 16 (PDF p.28): baseline TIN spectra vs Amiet predictions, same four Theta, three Phi,
    32 m/s, ~2x10^2 to 2x10^4 Hz, same SPL units (curves shifted by +5/-5/-15 dB as labeled).
- Shear-layer refraction: EXPLICITLY corrected. Amiet (1978) angle/amplitude abacuses
  (angle correction ~5 deg, amplitude 1-2 dB at extreme angles, Fig 12) plus a novel Kirchhoff
  integral method over the jet shear layers (Sec. 4.2). Fig 16 spectra use the Amiet abacuses.

**Vector or raster:** RASTER (one embedded PNG per figure, 1250-1460 px wide; good quality,
axis text legible at 150 dpi).

**Amiet comparison:** Fig 16 (PDF p.28): "A very good agreement is found, with discrepancies
below 1 dB over a wide range of frequencies and for most observation angles"; deviations only
below 300-500 Hz (mu ~ 0.3, near the high-frequency-theory validity limit, plus suspected
shear-layer oscillation contamination). Directivity: Figs 14-15 — corrected predictions in
"remarkably good agreement"; "Amiet's model is finally fully validated in a three-dimensional
context" (PDF p.27).

---

## Usability as ODT-RDT/Amiet validation cases (baseline spectra + U + Tu + Lambda + observer)

1. **Narayanan2015**: USABLE — baseline flat-plate SPL (Fig 4, 60 m/s) and PWL (Fig 5a, 40 m/s)
   with Tu 2.5%, Lambda 6 mm, R 1.2 m; caveat: Fig 4 observer angle NOT STATED; raster figures.
2. **Chaitanya2017**: MARGINAL — turbulence superbly documented (two length scales), but the only
   absolute baseline spectrum is NACA65 PWL vs fc0/U (Fig 10); no flat-plate absolute spectrum.
3. **Bampanis2019**: USABLE — baseline S_pp at 19/27/32 m/s, 90 deg, R 1.25 m, Tu 4.5%,
   Lambda 9 mm; vector figures; caveat: dB reference of S_pp axis not stated in this paper.
4. **Bampanis2022**: USABLE (best) — baseline TIN SPL in dB/Hz re 20 uPa at 3 speeds, multi-angle
   (4 azimuths x 3 elevations), stated shear-layer correction, and direct Amiet overlay (Fig 16).
5. Overall: Bampanis2022 is the primary quantitative anchor; Narayanan2015 and Bampanis2019 are
   secondary; Chaitanya2017 serves for serration/Delta-PWL physics, not absolute baselines.
