"""Step 2: ensemble eddy-on-mean injection rate, verified.

Claims tested (plane strain A=diag(a,-a,0), line along x2, U2 = -a*x2):

  C1  A single eddy of size l applied to the FULL velocity U2+u' changes the
      fluctuation energy on its interval by  DE = (4/27) a^2 l^3  on ensemble
      average, INDEPENDENT of the pre-existing turbulence level sigma
      (the cross term 2*int delta*M[u'] has zero mean), while an eddy applied
      to the fluctuation alone changes nothing at u'=0.
  C2  Domain-averaged injection over an event ensemble with size density p(l):
      dR22|ev = n_ev * (4/27) a^2 <l^3> / L   (bookkeeping/normalisation).
  C3  The discrete triplet map is an exact permutation (energy of a mapped
      field preserved to machine precision), and its action on the linear
      mean reproduces the continuum sawtooth variance (4/27)(a l)^2.

Consequence (derivation, PLAN_rebuild_eddymean.md): the eddy-on-mean channel
is quadratic in a and R-independent -- it does NOT reduce to the exact rapid
production P22 = 2 a R22 and therefore SUPPLEMENTS (with dose control), not
replaces, the continuous production.
"""
import numpy as np

rng = np.random.default_rng(7)

# ---------------- exact discrete triplet map (permutation) ----------------
def tm_perm(m3):
    """Permutation indices for the discrete triplet map on 0..m3-1 (m3%3==0):
    first third <- indices ==0 mod 3; middle third <- ==1 mod 3 reversed;
    last third <- ==2 mod 3."""
    i = np.arange(m3)
    return np.concatenate([i[0::3], i[1::3][::-1], i[2::3]])

def apply_eddy(v, i0, m3, perm=None):
    if perm is None:
        perm = tm_perm(m3)
    v = v.copy()
    v[i0:i0 + m3] = v[i0:i0 + m3][perm]
    return v

# ---------------- setup ----------------
L, N, a = 1.0, 30000, 1.0
dy = L / N
y = (np.arange(N) + 0.5) * dy
U = -a * y                                   # mean profile along the line

def smooth_noise(sigma, width=25):
    w = rng.standard_normal(N)
    k = np.exp(-0.5 * (np.arange(-3*width, 3*width+1) / width) ** 2)
    w = np.convolve(w, k / k.sum(), mode="same")
    s = w.std()
    return w * (sigma / s if s > 0 else 0.0)

print("=== C3: permutation exactness + sawtooth variance ===")
for m3 in (300, 3000, 9000):
    l = m3 * dy
    perm = tm_perm(m3)
    seg = U[:m3]
    d = seg[perm] - seg
    print(f"  l={l:5.3f}: var(delta)/(a l)^2 = {np.mean(d**2)/(a*l)**2:.6f} "
          f"(4/27={4/27:.6f});  sum(delta)*dy = {d.sum()*dy:+.2e}")
u_test = smooth_noise(1.0)
e0 = np.sum(u_test[:9000] ** 2)
e1 = np.sum(u_test[:9000][tm_perm(9000)] ** 2)
print(f"  energy preservation on random field: |dE|/E = {abs(e1-e0)/e0:.2e}")

print("\n=== C1: per-event injection, laminar and with turbulence ===")
m3 = 3000; l = m3 * dy; perm = tm_perm(m3)
pred = (4/27) * a**2 * l**3
print(f"  prediction DE = (4/27) a^2 l^3 = {pred:.4e}   (aloft: a*l = {a*l:.3f})")
for sigma_rel in (0.0, 0.5, 2.0, 5.0):
    sigma = sigma_rel * a * l
    K = 2500
    de = np.empty(K)
    for k in range(K):
        up = smooth_noise(sigma) if sigma > 0 else np.zeros(N)
        i0 = rng.integers(0, N - m3)
        v = U + up
        vn = v.copy(); vn[i0:i0+m3] = v[i0:i0+m3][perm]
        un = vn - U                       # new fluctuation
        de[k] = (np.sum(un[i0:i0+m3]**2) - np.sum(up[i0:i0+m3]**2)) * dy
    print(f"  sigma = {sigma_rel:3.1f}*(a l): <DE> = {de.mean():.4e} "
          f"+- {de.std(ddof=1)/np.sqrt(K):.1e}   ratio to pred = {de.mean()/pred:.3f}")
print("  fluctuation-only event at u'=0: DE = 0 identically (map of zero).")

print("\n=== C2: ensemble bookkeeping over a size distribution ===")
# sizes from p(l) ~ l^-2 on [l_a, l_b] (any density works for the identity)
l_a, l_b = 0.02, 0.30
n_ev = 4000
u_inv = rng.uniform(1/l_b, 1/l_a, n_ev)
ls = 1.0 / u_inv                              # p(l) ~ l^-2
m3s = np.maximum(3, (np.round(ls / dy / 3).astype(int) * 3))
tot = 0.0
for m3i in m3s:
    li = m3i * dy
    i0 = rng.integers(0, N - m3i)
    permi = tm_perm(m3i)
    seg = U[i0:i0+m3i]
    d = seg[permi] - seg
    tot += np.sum(d**2) * dy / L
pred_tot = (4/27) * a**2 * np.sum((m3s * dy)**3) / L
print(f"  measured  sum DR22 = {tot:.4e}")
print(f"  predicted sum (4/27) a^2 sum l^3 / L = {pred_tot:.4e}   "
      f"ratio = {tot/pred_tot:.4f}")
