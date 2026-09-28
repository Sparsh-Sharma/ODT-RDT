# Eddy-on-mean re-validation — T2 spectral verdict (2026-09-28)

Campaign: 6 × 1024-rlz fw ensembles on caro (eddymean-rebuild binary
2f7838a), gated (EM) + ungated-probe (EMU, LeddyMeanKt=1e9) + same-binary
baseline (I), at Sk/eps~0.8 (fw_S2) and ~16 (fw_S40); precursor protocol,
kvisc 1e-5. Analysis: `fw_eddymean.py` on the master pipeline (threeway.py
observables, 128^3 DNS + exact-RDT references, e=1, dump index 4).

## db22(k2) per band at e=1 (band edges k2(e)/kc(0) 0.30..12)

Sk/eps ~ 0.8:
  DNS   +nan +0.032 -0.033 -0.029 -0.032 -0.036 -0.040 -0.046   (target: bends NEGATIVE)
  I     +0.089 +0.103 +0.103 +0.096 +0.073 +0.075 +0.079 +0.076 (flat positive = known deficiency)
  EM    +0.105 +0.103 +0.066 +0.084 +0.073 +0.070 +0.067 +0.080 (= baseline; neutral)
  EMU   +0.147 +0.135 +0.120 +0.118 +0.100 +0.095 +0.086 +0.097 (MORE positive; wrong way)
  ODTK  +0.076 +0.094 +0.061 +0.054 +0.037 +0.040 +0.037 +0.024 (bends down toward DNS)

Sk/eps ~ 16:
  DNS   +nan +0.143 +0.007 -0.014 -0.037 -0.046 -0.049 -0.053
  I     +0.124 +0.114 +0.123 +0.115 +0.113 +0.111 +0.109 +0.107
  EM    +0.127 +0.117 +0.121 +0.120 +0.120 +0.115 +0.114 +0.116 (= baseline)
  EMU   +0.288 +0.254 +0.246 +0.251 +0.226 +0.203 +0.177 +0.158 (far MORE positive)
  ODTK  +0.123 +0.119 +0.115 +0.110 +0.108 +0.096 +0.097 +0.099

## Verdict: NEGATIVE (the spectral win is not there)

- DNS/RDT: db22 positive at large scales, crosses zero, NEGATIVE at fine
  scales (upwash below the transverse mean). That negative fine-scale
  branch is what the referee-1/3 "wavenumber-uniform" objection needs.
- Baseline SC-ODT: flat, positive, wavenumber-uniform (the deficiency).
- Gated eddy-on-mean: indistinguishable from baseline — no spectral effect
  (gate off once turbulent, by design). Consistent with the hs2 bands
  (all paired contrasts ~0).
- Ungated eddy-on-mean (ceiling): amplifies db22 uniformly MORE positive
  at every wavenumber (dramatically so at Sk/eps=16), moving AWAY from the
  DNS. It strengthens the energy-containing/upwash anisotropy uniformly; it
  does not manufacture the fine-scale negative branch.
- Only the archived adaptive-depth relaxation (ODTK) bends the fine scales
  down toward the DNS — it remains the paper's scale-local-relaxation win.

## Consequences
- T1 moments: neutral (established previously, gated recovers baseline).
- T2 spectral: NEGATIVE (this file).
- T3 interventions: eddy-on-mean cannot replace the adaptive-depth/clock
  relaxation — it moves db22 the wrong way, so it cannot sustain the fine
  scales in the DNS direction.
- Net: the rebuild's value is the TRANSITION-FROM-LAMINAR capability
  (Alan's foundational bar) at spectral neutrality (gated). It is NOT a
  spectral improvement and does not simplify the paper by dropping
  interventions. Fold in (if at all) as the foundational-correctness /
  transition result, not as a spectral advance. Decision + note-to-Alan
  belongs in the paper chat.
