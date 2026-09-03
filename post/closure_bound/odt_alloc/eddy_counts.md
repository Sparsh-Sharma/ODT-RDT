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

# Precursor campaign, S=2 cases (2026-09-03, tStrainOn = 0.4, job 4429881)

Accepted eddies per realization before onset and in each strain quarter
(t = 0.4 + e/2), mean post-onset eddy size (fraction of the line).

| case       | pre-onset | e∈[0,.25) | [.25,.5) | [.5,.75) | [.75,1] | total e≤1 | size |
|------------|----------:|----------:|---------:|---------:|--------:|----------:|-----:|
| S2_ISO     | 129.3 | 3.07 | 1.82 | 1.14 | 0.77 | 6.81 | 0.089 |
| S2_CHI     | 129.2 | 3.15 | 1.76 | 1.12 | 0.73 | 6.77 | 0.089 |
| S2_TYPESeq | 128.5 | 3.09 | 1.72 | 1.06 | 0.76 | 6.62 | 0.090 |
| S2_TYPESw  | 128.4 | 3.12 | 1.70 | 1.11 | 0.73 | 6.65 | 0.089 |

~6.7 relaxed eddies of size 0.09 L per line inside e≤1 (~0.6 line coverages),
as designed from the pilot; the burst (129 eddies) is entirely before onset.
Onset state at 1024 rlz: b = (-0.005, +0.003, +0.002) +- 0.003, i.e. isotropic.

# Precursor campaign, S=40 cases (rapid, Sk/eps = 16)

| case        | pre-onset | post-onset (t>=0.4) | max | realizations lost |
|-------------|----------:|--------------------:|----:|------------------:|
| S40_ISO     | 129.2 | 0.72 | 5 | 2 |
| S40_CHI     | 129.3 | 0.68 | 5 | 9 |
| S40_TYPESeq | 128.4 | 0.69 | 4 | 2 |
| S40_TYPESw  | 128.6 | 0.67 | 7 | 7 |

Under one eddy per line inside e<=1: a true rapid limit. Test 1 result: every
mode reproduces RDT (slope 0.129-0.133 vs 2/15, b22(1)-b22(0) = 0.120-0.124
vs 0.123), CHI included -- its S=2 deviation is the per-eddy beta bias.

Lost realizations (NaN-padded by the extractor, masked by good_rlz): they
abort in domain::domainPositionToIndex with an eddy position of about +-4 L,
far outside the line, only at S=40 (none at S=2). Robustness bug of the
line dilatation at very high strain rate; 0.2-0.9% of realizations, no
effect on the ensemble statistics. Open item.
