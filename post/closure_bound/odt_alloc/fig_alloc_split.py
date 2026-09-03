"""Figure for the allocation tests: (left) Test 1, the rapid-limit b_22(e)
per mode against the exact RDT value -- all modes must coincide; (right)
Test 3, the transverse splitting phi_11/phi_33 - 1 versus wavenumber at
S=0.5, e=1 per mode against the 128^3 DNS -- ISO predicted kappa-rigid,
CHI/TYPES kappa-dependent."""

import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from alloc_tests import (DNS_DIR, ECHK, MODES, SEG_FRAC, b_of_dump,  # noqa: E402
                         line_spectra, load)

COL = {"ISO": "#2a78d6", "CHI": "#eb6834", "TYPESeq": "#1baf7a",
       "TYPESw": "#eda100"}
MRK = {"ISO": "o", "CHI": "s", "TYPESeq": "^", "TYPESw": "D"}
LBL = {"ISO": "ISO (isotropic kernel, original)",
       "CHI": r"CHI ($\chi{=}0,\ \beta{=}0.3$)",
       "TYPESeq": "TYPES (equal $p_m$)",
       "TYPESw": "TYPES ($p_3{=}0.5$)"}
RDT_B22 = {0.25: 0.0383, 0.5: 0.0705, 1.0: 0.1229}

fig, (axl, axr) = plt.subplots(1, 2, figsize=(10.2, 3.8))

# ---- left: Test 1 ----
for mode in MODES:
    d = load(f"S8_{mode}")
    if d is None:
        continue
    b22, se = [], []
    for di in range(len(ECHK)):
        b, s = b_of_dump(d["lines"][di])
        b22.append(b[1])
        se.append(s[1])
    axl.errorbar(ECHK, b22, yerr=se, color=COL[mode], marker=MRK[mode],
                 ms=4.5, lw=1.5, capsize=2, label=LBL[mode])
axl.plot(list(RDT_B22), list(RDT_B22.values()), "k--", lw=1.2,
         label="RDT exact (Cauchy)")
axl.plot(ECHK, ECHK / 3.0, ":", color="0.5", lw=1.0,
         label="production only (slope 1/3)")
axl.set_xlabel(r"accumulated strain $e$", fontsize=9.5)
axl.set_ylabel(r"$b_{22}$", fontsize=9.5)
axl.set_title("Test 1: rapid limit ($S{=}8$) -- modes must coincide",
              fontsize=10)
axl.legend(fontsize=7.5, frameon=False, loc="upper left")

# ---- right: Test 3 ----
# Both systems on x = kappa_2(e)/kappa_p0, each normalized by its OWN initial
# peak (ODT: 8 waves per unit line; DNS: kp=4 in integer box units). Both
# dilate by exp(e/2), so the strain shift is directly comparable.
edges = np.geomspace(0.15, 8.0, 13)
dns = sorted(glob.glob(os.path.join(DNS_DIR, "chk_r0.8_s*_e1.npz")))
if dns:
    ph = np.mean([np.load(f, allow_pickle=True)["phi_line"] for f in dns],
                 axis=0)
    k2d = np.load(dns[0], allow_pickle=True)["kappa2"]
    use = (k2d > 0) & (k2d <= 0.85 * 42 * np.exp(0.5))
    xd = k2d / 4.0
    xc, yv = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (xd >= a) & (xd < b) & use
        if sel.sum():
            xc.append(np.sqrt(a * b))
            yv.append(ph[0][sel].sum() / ph[2][sel].sum() - 1)
    axr.plot(xc, yv, "k-", lw=2.2, label=r"DNS $128^3$ ($Sk/\varepsilon{=}0.8$)")
for mode in MODES:
    d = load(f"S05_{mode}")
    if d is None:
        continue
    k, (p11, _), (p33, _) = line_spectra(d["lines"][-1], d["Ldump"][-1])
    x = k / (2 * np.pi) / 8.0
    xc, yv = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (x >= a) & (x < b)
        if sel.sum():
            xc.append(np.sqrt(a * b))
            yv.append(p11[sel].sum() / p33[sel].sum() - 1)
    axr.plot(xc, yv, color=COL[mode], marker=MRK[mode], ms=4, lw=1.5,
             label=LBL[mode])
axr.axhline(0, color="0.75", lw=0.8, ls=":")
axr.set_xscale("log")
axr.set_xlabel(r"$\kappa_2(e)/\kappa_{p0}$", fontsize=9.5)
axr.set_ylabel(r"$\phi_{11}/\phi_{33}-1$", fontsize=9.5)
axr.set_title("Test 3: transverse splitting at $S{=}0.5$, $e{=}1$",
              fontsize=10)
axr.legend(fontsize=7.5, frameon=False, loc="lower left")

for ax in (axl, axr):
    ax.tick_params(labelsize=8.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.93", lw=0.6, zorder=0)
fig.suptitle("Kernel energy allocation: rapid response is mode-independent; "
             "the spectral splitting is where modes differ",
             fontsize=10.5, y=1.02)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, "fig_alloc_tests." + ext),
                bbox_inches="tight", dpi=180)
print("saved fig_alloc_tests.pdf/.png")
