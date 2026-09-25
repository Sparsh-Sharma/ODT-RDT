"""Two-panel figure for the note to Alan:
(a) an eddy of size l applied to the mean U2=-a x2 injects the sawtooth delta;
(b) the injected energy per event is independent of the turbulence level sigma
    (so eddy-on-mean supplements, not rescales, the fluctuation).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(7)
plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.8,
                     "figure.dpi": 130, "font.family": "DejaVu Sans"})

def tm_perm(m3):
    i = np.arange(m3)
    return np.concatenate([i[0::3], i[1::3][::-1], i[2::3]])

# ---------------- panel (a): the mechanism ----------------
L, N, a = 1.0, 30000, 1.0
dy = L / N
y = (np.arange(N) + 0.5) * dy
U = -a * (y - 0.5)                       # mean profile, centred for display
i0, m3 = int(0.35 * N), int(0.30 * N // 3 * 3)
l = m3 * dy
perm = tm_perm(m3)
MU = U.copy(); MU[i0:i0+m3] = U[i0:i0+m3][perm]
delta = MU - U

fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.4, 3.6))

axA.plot(y, U, color="0.45", lw=1.6, label=r"mean $U_2=-a\,x_2$")
axA.plot(y, MU, color="#1f5fa6", lw=1.6, label=r"$M[U_2]$ (after eddy)")
axA.plot(y, delta, color="#c0392b", lw=1.6, label=r"$\delta=(M-I)U_2$ (injected)")
axA.axvspan(i0*dy, (i0+m3)*dy, color="0.9", zorder=0)
axA.axhline(0, color="0.7", lw=0.6, zorder=0)
axA.annotate(f"eddy, size $\\ell$", xy=((i0+m3/2)*dy, -0.02),
             xytext=((i0+m3/2)*dy, -0.30), ha="center", fontsize=8.5,
             arrowprops=dict(arrowstyle="-", color="0.5", lw=0.7))
axA.set_xlabel(r"$x_2$ along the line"); axA.set_ylabel("velocity")
axA.set_title(r"(a) an eddy taps the mean", fontsize=10)
axA.legend(fontsize=8, loc="upper right", framealpha=0.9)
axA.text(0.03, 0.06, r"$\langle\delta^2\rangle=\frac{4}{27}(a\ell)^2$,  "
         r"$\langle\delta\rangle=0$", transform=axA.transAxes, fontsize=9.5,
         bbox=dict(boxstyle="round,pad=0.3", fc="#fff3f0", ec="#c0392b", lw=0.6))

# ---------------- panel (b): R-independence ----------------
def smooth_noise(sigma, width=25):
    w = rng.standard_normal(N)
    k = np.exp(-0.5 * (np.arange(-3*width, 3*width+1) / width) ** 2)
    w = np.convolve(w, k / k.sum(), mode="same")
    s = w.std()
    return w * (sigma / s if s > 0 else 0.0)

pred = (4/27) * a**2 * l**3
sig_rel = np.linspace(0, 6, 13)
K = 600
mean_r = np.empty_like(sig_rel); se_r = np.empty_like(sig_rel)
for j, sr in enumerate(sig_rel):
    sigma = sr * a * l
    de = np.empty(K)
    for k in range(K):
        up = smooth_noise(sigma) if sigma > 0 else np.zeros(N)
        ii = rng.integers(0, N - m3)
        un = up.copy(); un[ii:ii+m3] = up[ii:ii+m3][perm] + (U[ii:ii+m3][perm] - U[ii:ii+m3])
        de[k] = (np.sum(un[ii:ii+m3]**2) - np.sum(up[ii:ii+m3]**2)) * dy
    mean_r[j] = de.mean() / pred
    se_r[j] = de.std(ddof=1) / np.sqrt(K) / pred

axB.axhline(1.0, color="#c0392b", lw=1.2, ls="--",
            label=r"$\frac{4}{27}(a\ell)^2$ (turbulence-independent)")
axB.errorbar(sig_rel, mean_r, yerr=se_r, fmt="o", color="#1f5fa6", ms=4,
             capsize=2, lw=1, label="measured per-event injection")
axB.set_xlabel(r"pre-existing turbulence  $\sigma/(a\ell)$")
axB.set_ylabel(r"$\langle\Delta E\rangle\,/\,\frac{4}{27}(a\ell)^2$")
axB.set_title(r"(b) injection is independent of $R$", fontsize=10)
axB.set_ylim(0.6, 1.25); axB.legend(fontsize=8, loc="lower left")

fig.tight_layout()
fig.savefig("post/eddymean/fig_eddymean_note.png", dpi=130, bbox_inches="tight")
print("saved post/eddymean/fig_eddymean_note.png")
print("panel (b) ratios:", np.round(mean_r, 3))
