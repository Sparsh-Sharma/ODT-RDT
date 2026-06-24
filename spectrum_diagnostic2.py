#!/usr/bin/env python3
r"""
Spectrum diagnostic for the strain-coupled ODT (Level 1b, sec:val-spectrum).

For each dump it interpolates the (adapted, non-uniform, possibly dilated) line
onto a uniform grid, removes the mean, FFTs each velocity component, and forms
the 1-D component spectra E_i(k). It tracks two measures of how the spectrum
moves in k: the peak wavenumber and the spectral centroid kbar = <k E>/<E>.

The domain length is read FROM EACH DUMP, so the dilatation (line compression,
L(e) = L0 exp(A22 e)) is handled correctly: as L shrinks the wavenumbers
2*pi*n/L(e) grow, and the spectrum migrates to higher k.

VERIFICATION LADDER:
  (1) eddies OFF, dilatation OFF : spectrum only changes amplitude; peak and
      centroid stay put (ratio == 1). [already confirmed]
  (2) eddies OFF, dilatation ON  : pure linear RDT. NO cascade, but the geometric
      compression moves the peak by exactly exp(-A22 e). The dashed reference
      line is this prediction; the symbols must lie on it.
  (3) eddies ON,  dilatation ON  : genuine cascade adds transfer to high k, so
      the centroid rises ABOVE the linear-RDT line. That excess is the cascade.

Usage (Spyder): set DUMP_DIR, press Run.   CLI: python3 spectrum_diagnostic.py <dir>
"""
import os, glob, re, sys
import numpy as np
import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ======================================================================
DUMP_DIR = r'C:\Users\shar_sp\Documents\ODT-post\data\homogeneousStrain2'
SMAG     = 1.0          # strain magnitude S (plane a=0.5 -> S=1, e = t)
A22      = -0.5         # Astrain[2,2]; sets the linear-RDT peak migration exp(-A22 e)
NUNIFORM = 2048         # uniform resampling resolution for the FFT
STRAINS  = [0.0, 1.0, 2.0, 3.0, 4.0]
CENT_THRESH = 1e-3      # centroid integrates only where E>thresh*max (excludes interp floor)
OUT      = "fig_spectrum_evolution"
# ======================================================================

plt.rcParams.update({
    "text.usetex": False, "font.family": "serif",
    "font.serif": ["cmr10", "CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm", "axes.unicode_minus": False,
    "axes.formatter.use_mathtext": True,
    "font.size": 10, "axes.labelsize": 11, "axes.titlesize": 11,
    "legend.fontsize": 8.5, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "axes.linewidth": 0.8, "lines.linewidth": 1.4,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "figure.dpi": 150,
})
ccol = ["#1f77b4", "#d62728", "#2ca02c"]
cnm  = [r"$E_1(k)$", r"$E_2(k)$", r"$E_3(k)$"]

# ----------------------------------------------------------------------
def read_dump(fname):
    t=None; posf=[]; u=[]; v=[]; w=[]
    with open(fname) as f:
        for line in f:
            if line.startswith("#"):
                m=re.search(r"time\s*=\s*([-\d.eE+]+)", line)
                if m: t=float(m.group(1))
                continue
            if not line.strip(): continue
            c=line.split()
            posf.append(float(c[1])); u.append(float(c[2]))
            v.append(float(c[3])); w.append(float(c[4]))
    return t, np.array(posf), np.array(u), np.array(v), np.array(w)

def spectra_of(posf, fields, Nu):
    """Interpolate to a uniform grid over the actual (dilated) domain and FFT.
       Returns k[1:], [E1,E2,E3]. Domain is centred at 0: L = -2*posf[0]."""
    x0 = posf[0]; L = -2.0*x0
    faces = np.append(posf, -x0)                 # right boundary = -posf[0]
    xc = 0.5*(faces[:-1] + faces[1:])
    xu = x0 + (np.arange(Nu)+0.5)*(L/Nu)
    k  = 2.0*np.pi*np.fft.rfftfreq(Nu, d=L/Nu)
    out=[]
    for f in fields:
        fu = np.interp(xu, xc, f, period=L)
        fu = fu - fu.mean()
        fh = np.fft.rfft(fu)/Nu
        out.append(2.0*np.abs(fh)**2)
    return k[1:], [E[1:] for E in out]

def peak_and_centroid(k, E):
    pk = k[np.argmax(E)]
    m  = E > CENT_THRESH*E.max()                 # exclude interpolation floor
    kbar = (k[m]*E[m]).sum()/E[m].sum() if m.any() else np.nan
    return pk, kbar

# ----------------------------------------------------------------------
root = DUMP_DIR if len(sys.argv)<2 else sys.argv[1]
def find_dumps(d):
    for cand in (os.path.join(d,"data"), d):
        fs = glob.glob(os.path.join(cand,"data_*","dmp_*.dat")) or \
             glob.glob(os.path.join(cand,"dmp_*.dat"))
        if fs: return fs
    return []
files = find_dumps(root)
if not files:
    print(f"no dumps found under {root}"); sys.exit(1)

rows=[]
for f in files:
    t,posf,u,v,w = read_dump(f)
    if t is not None: rows.append((t,posf,u,v,w))
rows.sort(key=lambda r:r[0])
times = np.array([r[0]*SMAG for r in rows])
L0 = -2.0*rows[0][1][0]
print(f"read {len(rows)} dumps, e in [{times[0]:.2f}, {times[-1]:.2f}], L0 = {L0:.4f}")
print(f"final L/L0 = {(-2.0*rows[-1][1][0])/L0:.4f}  (exp(A22 e_max) = {np.exp(A22*times[-1]):.4f} if dilatation on)")

# ----------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.9))
peak={0:[],1:[],2:[]}; cent={0:[],1:[],2:[]}; epts=[]
for e_target in STRAINS:
    j=int(np.argmin(np.abs(times-e_target)))
    _,posf,u,v,w = rows[j]
    k,Es = spectra_of(posf,[u,v,w],NUNIFORM)
    shade=0.25+0.7*(e_target/max(STRAINS))
    for i in range(3):
        ax[0].loglog(k,Es[i],color=ccol[i],alpha=shade,lw=1.1)
    epts.append(times[j])
    for i in range(3):
        pk,kb = peak_and_centroid(k,Es[i]); peak[i].append(pk); cent[i].append(kb)
epts=np.array(epts)

ax[0].set_xlabel(r"wavenumber $k$"); ax[0].set_ylabel(r"$E_i(k)$")
ax[0].set_title(r"Component spectra (lighter$\to$darker $=$ increasing $e$)")
ax[0].grid(alpha=0.25,which="both",lw=0.4)
ax[0].legend([Line2D([0],[0],color=ccol[i],lw=2) for i in range(3)],cnm,
             loc="lower left",frameon=False)

# right: spectral centroid migration vs the linear-RDT prediction
ee=np.linspace(0,max(STRAINS),100)
ax[1].plot(ee,np.exp(-A22*ee),"k--",lw=1.2,label=r"linear RDT $e^{-A_{22}e}$")
for i in range(3):
    ax[1].plot(epts,np.array(cent[i])/cent[i][0],"o-",color=ccol[i],ms=3.6,
               label=fr"$\bar k_{i+1}/\bar k_{i+1}(0)$")
ax[1].axhline(1.0,color="k",lw=0.5,ls=":")
ax[1].set_xlabel(r"total strain $e=S\,t$")
ax[1].set_ylabel(r"spectral shift  $\bar k(e)/\bar k(0)$")
ax[1].set_title(r"Spectral migration vs linear RDT")
ax[1].grid(alpha=0.25,lw=0.4); ax[1].legend(frameon=False,fontsize=8)
fig.tight_layout()
for ext in ("pdf","png"): fig.savefig(f"{OUT}.{ext}",bbox_inches="tight")
print("wrote",os.path.abspath(f"{OUT}.pdf"),"/ .png")

print(f"\n{'e':>5}  peak_k1/0  cent_k1/0   RDT exp(-A22 e)")
for n,e in enumerate(epts):
    print(f"{e:5.2f}  {peak[0][n]/peak[0][0]:8.3f}  {cent[0][n]/cent[0][0]:8.3f}   {np.exp(-A22*e):8.3f}")
print("\n(2) eddies OFF + dilatation ON: centroid should track exp(-A22 e) (on the dashed line).")
print("(3) eddies ON  + dilatation ON: centroid rises ABOVE the dashed line -> cascade.")
plt.show()
