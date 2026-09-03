# Unstrained precursor pilot (2026-09-03, caro, 64 realizations, ISO, tStrainOn=1000)

Same IC and parameters as cases/odt_alloc/S05_ISO, never strained, 11 dumps.
k and eps on the native grid (eps = nu * sum_i <(du_i/dy)^2>, nu = 1e-4);
eddies from the runtime logs (accepted eddies per realization in [t, t+0.1),
mean size as a fraction of the line, coverage = rate * size).

| t    | k        | eps      | k/eps | b11    | b22    | b33    | eddies | rate  | size  | coverage/t |
|------|----------|----------|-------|--------|--------|--------|--------|-------|-------|------------|
| 0.00 | 9.96e-01 | 8.65e-01 | 1.15  | -0.001 | +0.000 | +0.001 | 96.8   | 968   | 0.032 | 30.9 |
| 0.10 | 2.59e-01 | 2.39e+00 | 0.11  | +0.004 | -0.009 | +0.005 | 19.5   | 195   | 0.049 | 9.6  |
| 0.20 | 1.20e-01 | 5.60e-01 | 0.21  | -0.000 | -0.011 | +0.011 | 8.1    | 81    | 0.062 | 5.0  |
| 0.30 | 7.55e-02 | 2.18e-01 | 0.35  | +0.003 | -0.016 | +0.013 | 4.7    | 47    | 0.073 | 3.4  |
| 0.40 | 5.36e-02 | 1.32e-01 | 0.41  | +0.002 | -0.018 | +0.016 | 2.4    | 24    | 0.075 | 1.8  |
| 0.50 | 4.15e-02 | 7.23e-02 | 0.58  | -0.001 | -0.023 | +0.024 | 2.1    | 21    | 0.086 | 1.8  |
| 0.60 | 3.39e-02 | 4.86e-02 | 0.70  | +0.001 | -0.024 | +0.023 | 1.8    | 18    | 0.098 | 1.7  |
| 0.70 | 2.83e-02 | 3.79e-02 | 0.75  | +0.002 | -0.026 | +0.024 | 1.1    | 11    | 0.107 | 1.2  |
| 0.80 | 2.42e-02 | 2.71e-02 | 0.89  | +0.003 | -0.027 | +0.024 | 1.2    | 12    | 0.109 | 1.3  |
| 0.90 | 2.10e-02 | 2.05e-02 | 1.02  | +0.008 | -0.028 | +0.020 | 0.5    | 5     | 0.115 | 0.6  |
| 1.00 | 1.86e-02 | 1.40e-02 | 1.33  | +0.010 | -0.029 | +0.019 |        |       |       |      |

Reading: the IC-relaxation burst (97 eddies, 75% of k lost) is over by
t~0.2; afterwards k/eps ~ 1.2 t and the eddy coverage time 1/(rate*size)
agrees with k/eps to 10-30%, so the ODT analogue of Sk/eps is unambiguous.
The unstrained line drifts to b22 ~ -0.03, b33 ~ +0.02 (64 rlz; to be
confirmed at 1024 from the onset dump) -- tests must use b_ij(e) - b_ij(0).

Campaign choice: tStrainOn = 0.4 (k/eps = 0.41, coverage time 0.55);
S = 2  -> S k/eps = 0.8  (DNS slow case),  e = 1 at t = 0.9,   ~8.5 eddies/rlz in e<=1 (~0.8 line coverages);
S = 40 -> S k/eps = 16   (DNS rapid case), e = 1 at t = 0.425, ~0.6 eddies/rlz.
Dumps at e = 0, .25, .5, .75, 1 (explicit list).
