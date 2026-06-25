# Section 4.2 — Emergent distorted spectrum: clean run protocol

One case folder, four runs, each a single switch away from the last. Do them in
order; each has a PASS check so you know it worked before moving on. Everything
uses the SAME base input (`input/strainSpectrum/input.yaml`); each run only flips
the flags listed. Use a fresh case name per run so nothing overwrites.

Base physical setup (do NOT change between runs unless told):
  probType      HOMOGENEOUS_STRAIN
  domainLength  1.0          (non-dimensional; nu = kvisc0 ~ 1/Re)
  Astrain       diag(0.5,-0.5,0)   -> S=1, e=t ; A22=-0.5 ; line compresses as L=L0 e^(A22 e)
  tEnd          4.0
  ngrd0         4000
  strainClosure LRR
  spectral IC   specKpWaves 8 , specNmodes 64
  strainCFL     0.01         (dt set by strain, NOT dump cadence)
  Lstrain       true         (always on for all four runs)
  dumpTimesGen  0 -> 4 step 0.1   (40 dumps)

====================================================================
RUN A1 — rapid-distortion limit (moments check)        case: sp_A1_rdt
====================================================================
Purpose: confirm the implementation still reproduces the LRR moment trajectory
         (this is the Section 4.1.5 check, repeated as the spectrum baseline).
Flags:   LnoEddies   true
         Ldilatation false
         kvisc0      0.0
Script:  analyze_sweep.py   (point at sp_A1_rdt)
PASS:    fractions on LRR (0.097 / 0.601 / 0.302), kt deficit ~ -0.3%, 40 dumps.
Figure:  (already have it as fig:level1a — no new figure needed here.)

====================================================================
RUN A2 — pure compression (linear-RDT spectral test)   case: sp_A2_dil
====================================================================
Purpose: the dilatation alone must reproduce linear-RDT spectral compression:
         peak/centroid migrate by exactly exp(-A22 e), NO cascade, moments
         unchanged from A1.
Flags:   LnoEddies   true
         Ldilatation true
         kvisc0      0.0
Scripts: spectrum_diagnostic.py (point at sp_A2_dil, single-case mode)
         analyze_sweep.py       (confirm moments still on LRR)
PASS:    centroid lands ON the dashed exp(-A22 e) line -> 7.39 at e=4
         (you saw 7.029 at e=3.9, exact); final L/L0 = exp(A22*e_max) ~ 0.142;
         fractions still 0.097/0.601/0.302.
Figure:  fig:spec-rdt  (single-case spectrum_diagnostic output)

====================================================================
RUN B — full model, eddies + compression               case: sp_B_full
====================================================================
Purpose: the emergent distorted spectrum — eddies add a cascade ON TOP of the
         compression. THIS NEEDS THE BOUNDARY-CLAMP FIX in domain.cc
         (domainPositionToIndex) or it crashes at e~0.3.
Flags:   LnoEddies   false
         Ldilatation true
         kvisc0      <pick physically, see RUN C sweep; start 1.0e-5>
PASS:    runs to e=4 (40 dumps); centroid sits ABOVE the A2 compression line at
         the resolved viscosities; spectrum broadens beyond the A2 shape.
Figure:  fig:spec-cascade (single-case spectrum_diagnostic output at chosen nu)

====================================================================
RUN C — Reynolds-number (kvisc0) sweep                  cases: hs2_kv*
====================================================================
Purpose: map cascade depth vs nu and locate the physical+resolved window.
Flags:   LnoEddies   false ; Ldilatation true  (forced by the script per case)
         kvisc0      swept: 1e-4, 5e-5, 2e-5, 1.42e-5, 1e-5, 1e-6, 1e-7
Scripts: run_kvisc_sweep.sh   (parallel; sets flags per case)
         analyze_kvisc_sweep.py   (centroid-vs-RDT table)
         spectrum_diagnostic.py   (point at PARENT data/ -> comparison figure)
PASS/READ:
   - centroid climbs monotonically as nu drops (table).
   - BUT low-nu (1e-6, 1e-7) show a FLAT high-k shelf at the Nyquist =
     grid-scale pileup (under-resolved); their centroid is partly numerical.
   - resolved cases (1e-5) sit BELOW the RDT line; no clean k^-5/3 yet.
   - CONCLUSION: cascade deepens monotonically with Re (correct physics, beyond
     linear RDT), but quantitative depth is resolution-limited at low nu.
Figures: fig:spec-compare (comparison spectrum_diagnostic over all kv cases)

====================================================================
RESOLUTION CHECK (do before quoting any low-nu number)  case: sp_B_hires
====================================================================
Take nu=1e-6, rerun with ngrd0=16000, dxmin=5e-5 (all else as RUN B).
If the high-k shelf drops / rolloff steepens -> pileup confirmed, that nu needs
this resolution. If unchanged -> cascade is genuine.

====================================================================
HOW TO PICK kvisc0 PHYSICALLY (not from the centroid)
====================================================================
nu = kvisc0 ~ 1/Re_L of the TARGET inflow turbulence (experiment or LES/LBM you
match for Gate 2). 1e-4 -> Re_L~1e4 ; 1e-6 -> Re_L~1e6. Choose nu = 1/Re_target,
THEN verify with spectrum_diagnostic (single mode) that the rolloff sits below
the Nyquist with margin. Do NOT pick nu to make the centroid beat RDT.

====================================================================
DIRECTORY HYGIENE (end the hotch-potch)
====================================================================
data/sp_A1_rdt , data/sp_A2_dil , data/sp_B_full , data/sp_B_hires , data/hs2_kv*
Delete everything else (homogeneousStrain, *1b, *1c, *2, stray hs_*). Each run =
its own caseName; never reuse a folder.
