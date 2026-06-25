#!/usr/bin/env python3
r"""
Spectrum diagnostic for the strain-coupled ODT (Level 1b, sec:val-spectrum).

TWO MODES, auto-detected from DUMP_DIR:

  * SINGLE CASE  -- DUMP_DIR points at one case folder (has data/data_*/dmp_*.dat):
        left  : component spectra E_i(k) at several strains (lighter->darker)
        right : spectral centroid migration vs the linear-RDT line exp(-A22 e)

  * COMPARISON   -- DUMP_DIR points at a parent holding several case folders
                    (e.g. data/ with hs2_kv1em4, hs2_kv1em5, ...):
        left  : final-strain total spectrum E(k)=Sum_i E_i for each case, overlaid,
                with a -5/3 reference slope -> shows the inertial range developing
                as kvisc0 drops (the cascade getting cleaner / deeper)
        right : centroid migration cbar(e_max)/cbar(0) per case vs the RDT line

The domain length is read FROM EACH DUMP, so the dilatation (L(e)=L0 exp(A22 e))
is handled: as L shrinks the wavenumbers 2*pi*n/L grow and the spectrum migrates.

Usage (Spyder): set DUMP_DIR, press Run.   CLI: python3 spectrum_diagnostic.py <dir>
"""
import os, glob, re, sys
import numpy as np
import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ======================================================================
# Point at ONE case for detail, or at the PARENT data/ dir for comparison:
DUMP_DIR = r'\\wsl.localhost\Ubuntu\home\shar_sp\ODT-RDT\data'
SMAG     = 1.0          # strain magnitude S (plane a=0.5 -> S=1, e = t)
A22      = -0.5         # Astrain[2,2]; linear-RDT migration is exp(-A22 e)
NUNIFORM = 4096         # uniform resampling resolution for the FFT
STRAINS  = [0.0, 1.0, 2.0, 3.0, 4.0]    # single-case mode: strains to show
CENT_THRESH = 1e-3      # centroid integrates only where E>thresh*max (drops floor)
OUT      = "fig_spectrum"
# ======================================================================

plt.rcParams.update({
    "text.usetex": False, "font.family": "serif",
    "font.serif": ["cmr10", "CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm", "axes.unicode_minus": False,
    "axes.formatter.use_mathtext": True,
    "font.size": 10, "axes.labelsize": 11, "axes.titlesize": 11,
    "legend.fontsize": 8, "xtick.labelsize": 9, "ytick.labelsize": 9,
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
    rows.sort(key=lambda r:r[0])
    return rows

def nice_label(name):
    # hs2_kv1em6 -> "kvisc0=1e-6"; hs2_kv1_42em5 -> "1.42e-5"
    s=name.replace("hs2_","").replace("kv","")
    s=s.replace("_",".").replace("em","e-").replace("ep","e+")
    return f"$\\nu$={s}"

# ----------------------------------------------------------------------
root = DUMP_DIR if len(sys.argv)<2 else sys.argv[1]

# discover cases: is root itself a case, or a parent of cases?
if dumps_in(root):
    cases=[root]
else:
    cases=sorted(d for d in glob.glob(os.path.join(root,"*"))
                 if os.path.isdir(d) and dumps_in(d))
if not cases:
    print(f"no dumps found at/under {root}"); sys.exit(1)

# ======================================================================
if len(cases)==1:
    # -------- SINGLE-CASE DETAIL --------
    rows=load_case(cases[0]); times=np.array([r[0]*SMAG for r in rows])
    L0=-2.0*rows[0][1][0]
    print(f"[single] {os.path.basename(cases[0])}: {len(rows)} dumps, "
          f"e in [{times[0]:.2f},{times[-1]:.2f}], final L/L0={(-2*rows[-1][1][0])/L0:.4f}")
    fig,ax=plt.subplots(1,2,figsize=(9.4,3.9))
    peak={0:[],1:[],2:[]}; cent={0:[],1:[],2:[]}; epts=[]
    for e_t in STRAINS:
        j=int(np.argmin(np.abs(times-e_t))); _,posf,u,v,w=rows[j]
        k,Es=spectra_of(posf,[u,v,w],NUNIFORM); shade=0.25+0.7*(e_t/max(STRAINS))
        for i in range(3): ax[0].loglog(k,Es[i],color=ccol[i],alpha=shade,lw=1.1)
        epts.append(times[j])
        for i in range(3):
            cent[i].append(centroid(k,Es[i])); peak[i].append(k[np.argmax(Es[i])])
    epts=np.array(epts)
    ax[0].set_xlabel(r"wavenumber $k$"); ax[0].set_ylabel(r"$E_i(k)$")
    ax[0].set_title(r"Component spectra (lighter$\to$darker $=$ increasing $e$)")
    ax[0].grid(alpha=0.25,which="both",lw=0.4)
    ax[0].legend([Line2D([0],[0],color=ccol[i],lw=2) for i in range(3)],cnm,
                 loc="lower left",frameon=False)
    ee=np.linspace(0,max(STRAINS),100)
    ax[1].plot(ee,np.exp(-A22*ee),"k--",lw=1.2,label=r"linear RDT $e^{-A_{22}e}$")
    for i in range(3):
        ax[1].plot(epts,np.array(cent[i])/cent[i][0],"o-",color=ccol[i],ms=3.6,
                   label=fr"$\bar k_{i+1}/\bar k_{i+1}(0)$")
    ax[1].set_xlabel(r"total strain $e=S\,t$"); ax[1].set_ylabel(r"$\bar k(e)/\bar k(0)$")
    ax[1].set_title(r"Spectral migration vs linear RDT")
    ax[1].grid(alpha=0.25,lw=0.4); ax[1].legend(frameon=False,fontsize=8)
    fig.tight_layout()
    for ext in ("pdf","png"): fig.savefig(f"{OUT}_single.{ext}",bbox_inches="tight")
    print("wrote",os.path.abspath(f"{OUT}_single.pdf"))

else:
    # -------- MULTI-CASE COMPARISON --------
    print(f"[compare] {len(cases)} cases under {root}")
    cmap=plt.cm.viridis(np.linspace(0,0.92,len(cases)))
    fig,ax=plt.subplots(1,2,figsize=(9.8,3.9))
    for c,col in zip(cases,cmap):
        rows=load_case(c)
        if len(rows)<2: print(f"  skip {os.path.basename(c)} (<2 dumps)"); continue
        times=np.array([r[0]*SMAG for r in rows])
        # final-strain total spectrum
        _,pN,uN,vN,wN=rows[-1]; k,Es=spectra_of(pN,[uN,vN,wN],NUNIFORM)
        Etot=Es[0]+Es[1]+Es[2]
        ax[0].loglog(k,Etot,color=col,lw=1.1,label=nice_label(os.path.basename(c)))
        # centroid migration (total spectrum) per case
        cN=[]; ee=[]
        for (t,p,u,v,w) in rows:
            kk,E=spectra_of(p,[u,v,w],NUNIFORM); cN.append(centroid(kk,E[0]+E[1]+E[2])); ee.append(t*SMAG)
        cN=np.array(cN); ee=np.array(ee)
        ax[1].plot(ee,cN/cN[0],"-",color=col,lw=1.3,label=nice_label(os.path.basename(c)))
    # -5/3 reference slope on the spectra panel
    kref=np.array([3e1,3e2]); ax[0].loglog(kref,2e-1*(kref/kref[0])**(-5/3),
                 "k:",lw=1.0); ax[0].text(1e2,3e-2,r"$k^{-5/3}$",fontsize=9)
    ax[0].set_xlabel(r"wavenumber $k$"); ax[0].set_ylabel(r"$E(k)=\sum_i E_i(k)$")
    ax[0].set_title(r"Final-strain spectrum vs $\nu$ (cascade depth)")
    ax[0].grid(alpha=0.25,which="both",lw=0.4); ax[0].legend(frameon=False,fontsize=7.5,loc="lower left")
    ee=np.linspace(0,max(STRAINS),100)
    ax[1].plot(ee,np.exp(-A22*ee),"k--",lw=1.3,label=r"linear RDT")
    ax[1].set_xlabel(r"total strain $e=S\,t$"); ax[1].set_ylabel(r"$\bar k(e)/\bar k(0)$")
    ax[1].set_title(r"Centroid migration vs $\nu$"); ax[1].grid(alpha=0.25,lw=0.4)
    ax[1].legend(frameon=False,fontsize=7.5,loc="upper left")
    fig.tight_layout()
    for ext in ("pdf","png"): fig.savefig(f"{OUT}_compare.{ext}",bbox_inches="tight")
    print("wrote",os.path.abspath(f"{OUT}_compare.pdf"))

plt.show()