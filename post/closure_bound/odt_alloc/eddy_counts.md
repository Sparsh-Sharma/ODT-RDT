# Accepted-eddy statistics of the allocation ensembles (2026-09-03)

Source: `data/<case>/runtime/runtime_*` on caro (one row per accepted eddy),
1024 realizations per case, strain from t=0, no precursor. Columns: mean
eddies per realization before t=0.05, cumulative mean count at each strain
checkpoint e (t = e/S), mean eddy size (fraction of the line) overall and
after t=0.1.

| case        | t<0.05 | e≤0.25 | e≤0.5 | e≤0.75 | e≤1.0 | size | size(t>0.1) |
|-------------|-------:|-------:|------:|-------:|------:|-----:|------------:|
| S8_ISO      |  64.4  |  46.5  | 72.3  |  83.8  |  89.4 | 0.031 | 0.041 |
| S8_CHI      |  63.9  |  46.3  | 71.6  |  82.9  |  88.7 | 0.031 | 0.040 |
| S8_TYPESeq  |  64.7  |  46.9  | 72.4  |  83.7  |  89.1 | 0.031 | 0.040 |
| S8_TYPESw   |  64.8  |  47.0  | 72.5  |  83.7  |  89.2 | 0.031 | 0.040 |
| S05_ISO     |  68.3  | 128.5  | 133.3 | 134.7  | 135.4 | 0.042 | 0.066 |
| S05_CHI     |  68.2  | 129.4  | 134.3 | 135.7  | 136.4 | 0.042 | 0.066 |
| S05_TYPESeq |  68.8  | 127.8  | 132.4 | 133.8  | 134.5 | 0.041 | 0.065 |
| S05_TYPESw  |  68.9  | 127.7  | 132.4 | 133.9  | 134.6 | 0.041 | 0.066 |

## Reading

* Every case starts with the same burst: ~65 eddies before t=0.05 (initial
  condition relaxing to an ODT-consistent state), independent of S and mode.
* S=8: all 89 eddies inside e≤1 belong to that burst (e=1 is t=0.125). The
  "rapid" runs are burst + continuous rapid operator, not a clean rapid limit.
* S=0.5: 128 eddies by e=0.25 (t=0.5), then only ~7 more over e in [0.25,1]
  (t from 0.5 to 2). The line is eddy-quiescent during the strained
  evolution, so the slow term is effectively absent in the "slow" runs too.
* Consequences: ISO b22(e) is the same curve at S=8 and S=0.5 (0.103 vs
  0.107 at e=1, jackknife SE 0.002/0.004); the transverse splitting is
  kappa-rigid for every mode (continuous-B response of a quiescent line);
  the only mode that separates is CHI, whose strain-aligned beta term is a
  per-eddy static bias acting on the burst, not a slow-term effect.
* Post-burst eddy rate at S=0.5 is ~5 per line per unit time with eddy size
  ~0.066 L, i.e. an eddy turnover ~3 time units, so S*tau_eddy ~ 1.5 versus
  the DNS Sk/eps = 0.8 and, more importantly, ~0.7 turnovers within e=1
  versus ~2.5 in the DNS.

## What a discriminating test needs

1. An unstrained ODT precursor long enough for the burst to pass and the
   eddy rate to settle (t ~ 0.2-0.3 on this IC), with the strain switched on
   afterwards (requires a strain-onset time in the code or a restart from
   the relaxed dump).
2. S chosen from the measured post-burst eddy rate so that S*tau_eddy
   matches the DNS Sk/eps = 0.8 (roughly S ~ 0.25 on this IC, to be
   re-measured after the precursor), giving O(2-3) eddy turnovers within e=1.
