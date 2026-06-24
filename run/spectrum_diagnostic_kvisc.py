# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 19:00:54 2026

@author: shar_sp
"""

#!/usr/bin/env python3
r"""
Comparison spectrum diagnostic (Section 4.2) for a kvisc0 sweep.

Point DUMP_DIR at the PARENT data/ directory holding the case folders
(hs2_kv1em4, hs2_kv1em5, ...). Produces TWO separate, manuscript-ready
figures, each as .pdf AND .png:

   fig_spectrum_compare_spectra   final-strain spectrum E(k)=Sum_i E_i vs nu
   fig_spectrum_compare_centroid  spectral-centroid migration per case vs nu

CENTROID SMOOTHING
   Each centroid point is from a single ODT line (one realization), so with
   eddies on the curve is jagged from stochastic eddy sampling -- that scatter
   is real, not a plotting artefact. SMOOTH_WIN applies a centred moving average
   of that many dumps in strain. Set SMOOTH_WIN=1 to see the raw data; if you
   smooth, SAY SO in the caption. The physically correct fix is to ENSEMBLE
   AVERAGE several seeds per case (see the multi-seed note at the bottom).
"""
import os, glob, re, sys
import numpy as np
import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ======================================================================
DUMP_DIR    = r'\\wsl.localhost\Ubuntu\home\shar_sp\ODT-RDT\data'
SMAG        = 1.0
A22         = -0.5
NUNIFORM    = 4096
CENT_THRESH = 1e-3
SMOOTH_WIN  = 5         # centred moving-average window (dumps). 1 = raw, no smoothing.
OUTBASE     = "fig_spectrum_compare"
SAVE_EXT    = ("pdf", "png")
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

def savefig(fig, name):
    for ext in SAVE_EXT: fig.savefig(f"{name}.{ext}", bbox_inches="tight")
    print("wrote", os.path.abspath(f"{name}.pdf"), "(+ .png)")

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
            posf.append(float(c[1])); u.append(float(c[2])); v.append(float(c[3])); w.append(float(c[4]))
    return t, np.array(posf), np.array(u), np.array(v), np.array(w)

def spectra_of(posf, fields, Nu):
    x0=posf[0]; L=-2.0*x0
    faces=np.append(posf,-x0); xc=0.5*(faces[:-1]+faces[1:])
    xu=x0+(np.arange(Nu)+0.5)*(L/Nu)
    k=2.0*np.pi*np.fft.rfftfreq(Nu,d=L/Nu)
    out=[]
    for f in fields:
        fu=np.interp(xu,xc,f,period=L); fu=fu-fu.mean()
        out.append(2.0*np.abs(np.fft.rfft(fu)/Nu)**2)
    return k[1:], [E[1:] for E in out]

def centroid(k, E):
    m=E>CENT_THRESH*E.max()
    return (k[m]*E[m]).sum()/E[m].sum() if m.any() else np.nan

def smooth(y, win):
    if win<=1 or y.size<3: return y
    win=min(win, y.size if y.size%2 else y.size-1)
    if win%2==0: win-=1
    if win<3: return y
    pad=win//2
    yp=np.pad(y, pad, mode="edge")
    ker=np.ones(win)/win
    return np.convolve(yp, ker, mode="valid")

def dumps_in(case):
    for cand in (os.path.join(case,"data"), case):
        fs=glob.glob(os.path.join(cand,"data_*","dmp_*.dat")) or \
           glob.glob(os.path.join(cand,"dmp_*.dat"))
        if fs: return fs
    return []

def load_case(case):
    rows=[]
    for f in dumps_in(case):
        t,posf,u,v,w=read_dump(f)
        if t is not None: rows.append((t,posf,u,v,w))
    rows.sort(key=lambda r:r[0]); return rows

def nice_label(name):
    s=name.replace("hs2_","").replace("kv","").replace("_",".").replace("em","e-").replace("ep","e+")
    return f"$\\nu$={s}"

# ----------------------------------------------------------------------
root = DUMP_DIR if len(sys.argv)<2 else sys.argv[1]
if dumps_in(root):
    cases=[root]
else:
    cases=sorted(d for d in glob.glob(os.path.join(root,"*"))
                 if os.path.isdir(d) and dumps_in(d))
if not cases:
    print(f"no dumps found at/under {root}"); sys.exit(1)
print(f"[compare] {len(cases)} cases under {root}  (SMOOTH_WIN={SMOOTH_WIN})")
cmap=plt.cm.viridis(np.linspace(0,0.92,len(cases)))

# ---- FIGURE 1: final-strain spectra vs nu ----
fig1,ax=plt.subplots(figsize=(5.8,4.2))
centro={}
for c,col in zip(cases,cmap):
    rows=load_case(c)
    if len(rows)<2: print(f"  skip {os.path.basename(c)} (<2 dumps)"); continue
    _,pN,uN,vN,wN=rows[-1]; k,Es=spectra_of(pN,[uN,vN,wN],NUNIFORM)
    ax.loglog(k,Es[0]+Es[1]+Es[2],color=col,lw=1.1,label=nice_label(os.path.basename(c)))
    cN=[]; ee=[]
    for (t,p,u,v,w) in rows:
        kk,E=spectra_of(p,[u,v,w],NUNIFORM); cN.append(centroid(kk,E[0]+E[1]+E[2])); ee.append(t*SMAG)
    cN=np.array(cN); ee=np.array(ee)
    centro[os.path.basename(c)]=(ee, smooth(cN,SMOOTH_WIN)/cN[0], col)
kref=np.array([3e1,3e2]); ax.loglog(kref,2e-1*(kref/kref[0])**(-5/3),"k:",lw=1.0)
ax.text(1e2,3e-2,r"$k^{-5/3}$",fontsize=9)
ax.set_xlabel(r"wavenumber $k$"); ax.set_ylabel(r"$E(k)=\sum_i E_i(k)$")
ax.set_title(r"Final-strain spectrum vs $\nu$ (cascade depth)")
ax.grid(alpha=0.25,which="both",lw=0.4); ax.legend(frameon=False,fontsize=8,loc="lower left")
fig1.tight_layout(); savefig(fig1,f"{OUTBASE}_spectra")

# ---- FIGURE 2: centroid migration vs nu ----
fig2,ax=plt.subplots(figsize=(5.8,4.2))
ee=np.linspace(0,4.0,100)
ax.plot(ee,np.exp(-A22*ee),"k--",lw=1.4,label=r"linear RDT")
for name,(eex,cc,col) in centro.items():
    ax.plot(eex,cc,"-",color=col,lw=1.5,label=nice_label(name))
ax.set_xlabel(r"total strain $e=S\,t$"); ax.set_ylabel(r"$\bar k(e)/\bar k(0)$")
ttl=r"Centroid migration vs $\nu$"
if SMOOTH_WIN>1: ttl += fr" ({SMOOTH_WIN}-pt moving avg.)"
ax.set_title(ttl)
ax.grid(alpha=0.25,lw=0.4); ax.legend(frameon=False,fontsize=8,loc="upper left")
fig2.tight_layout(); savefig(fig2,f"{OUTBASE}_centroid")

plt.show()

# ----------------------------------------------------------------------
# MULTI-SEED (ensemble) note: for the final figure, run each case a few times
# with different `seed` in the input (e.g. seeds 22,23,24,25), put them under
# data/<case>/seed_XX/, average cN over seeds at each strain, and plot the mean
# (optionally with a +/- std band). That removes the scatter physically rather
# than by smoothing, and is the version to put in the manuscript.