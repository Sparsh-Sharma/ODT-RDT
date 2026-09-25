"""Step 1 foundation: does an ODT triplet map on the mean profile U2 = -a*x2
inject a fluctuation of energy ~ (a*l)^2 at the eddy scale l?

Standard measure-preserving, continuous triplet map on [y0, y0+l], in
normalized coordinate s=(y-y0)/l with pre-image sigma(s):
  s in [0,1/3]  : sigma = 3s
  s in [1/3,2/3]: sigma = 2 - 3s   (reversed middle copy -> continuity)
  s in [2/3,1]  : sigma = 3s - 2
Verifies against the analytic result  <delta^2> = (4/27) (a l)^2 , <delta>=0.
"""
import numpy as np

def preimage_sigma(s):
    sig = np.empty_like(s)
    m1 = s < 1/3
    m2 = (s >= 1/3) & (s < 2/3)
    m3 = s >= 2/3
    sig[m1] = 3*s[m1]
    sig[m2] = 2 - 3*s[m2]
    sig[m3] = 3*s[m3] - 2
    return sig

def triplet_map(f_of_y, y0, l, N=30000):
    """Apply the triplet map to a field f (callable of y) on [y0,y0+l];
    return y grid, mapped field, original field on the interval."""
    s = (np.arange(N) + 0.5) / N          # cell centres in [0,1]
    y = y0 + s*l
    sigma = preimage_sigma(s)
    y_pre = y0 + sigma*l
    return y, f_of_y(y_pre), f_of_y(y)

def check(a, l, y0=0.0):
    f = lambda y: -a*y                    # mean profile U2 = -a x2
    y, Mf, f0 = triplet_map(f, y0, l)
    delta = Mf - f0                       # injected fluctuation
    ds = 1.0/len(y)
    mean_delta = delta.mean()
    var_delta = (delta**2).mean()         # = (1/l) int delta^2 dy  (s-average)
    # measure preservation: integral of Mf over interval == integral of f0
    int_Mf = Mf.mean()*l; int_f0 = f0.mean()*l
    return mean_delta, var_delta, int_Mf, int_f0

print("=== single eddy, a=1, l=1 ===")
md, var, iMf, if0 = check(1.0, 1.0)
print(f"  <delta>            = {md:+.3e}   (expect 0)")
print(f"  <delta^2>/(a l)^2  = {var:.6f}   (expect 4/27 = {4/27:.6f})")
print(f"  RMS/(a l)          = {np.sqrt(var):.6f}   (expect {np.sqrt(4/27):.6f})")
print(f"  measure preserved  = int Mf {iMf:+.4e} vs int f {if0:+.4e}")

print("\n=== energy scaling with eddy size l (a=1) ===")
print("   l        <delta^2>      /(a l)^2")
for l in [0.25, 0.5, 1.0, 2.0, 4.0]:
    _, var, _, _ = check(1.0, l)
    print(f"  {l:5.2f}   {var:.6e}   {var/(1.0*l)**2:.6f}")

print("\n=== linearity in strain rate a (l=1) ===")
print("   a        RMS(delta)     /(a l)")
for a in [0.5, 1.0, 2.0]:
    _, var, _, _ = check(a, 1.0)
    print(f"  {a:5.2f}   {np.sqrt(var):.6e}   {np.sqrt(var)/(a*1.0):.6f}")

# spectral content: where in wavenumber does the injected energy sit?
print("\n=== scale content of the injected fluctuation (l=1) ===")
y, Mf, f0 = triplet_map(lambda y: -1.0*y, 0.0, 1.0, N=3000)
d = Mf - f0
d = d - d.mean()
F = np.fft.rfft(d)/len(d)
k = np.fft.rfftfreq(len(d), d=1.0/len(d))   # in cycles per interval length l
P = np.abs(F)**2
kpk = k[1:][np.argmax(P[1:])]
print(f"  dominant wavenumber = {kpk:.1f} cycles per l  (fundamental of the sawtooth ~ 1.5/l)")
print(f"  fraction of energy in k<=3 (large scale ~l): {P[1:4].sum()/P[1:].sum():.3f}")
