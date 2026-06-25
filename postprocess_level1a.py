#!/usr/bin/env python3
r"""
Post-process a strain-coupled ODT (Level 1a) run: read the dmp_*.dat series,
compute the length-weighted component-energy fractions u_i^2/2k_t at each dump
time, and plot them against total strain e = S t, overlaid on the Level 0
benchmark (exact RDT and the LRR closure). Output matches the paper figure
style (Computer Modern, colour = component).

Usage:
    python3 postprocess_level1a.py [DUMP_DIR]
DUMP_DIR defaults to the current directory; it is scanned for dmp_*.dat
(and odt_init.dat / odt_end.dat if present). Strain rate is normalised so
S = 1 (plane strain a = 1/2), hence e = t; adjust SMAG if you change Astrain.
"""
import sys, glob, os, re
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

DUMP_DIR = sys.argv[1] if len(sys.argv) > 1 else "."
SMAG     = 1.0          # strain magnitude S used in the run (plane a=0.5 -> S=1)
OUT      = "fig_level1a_trajectory"

# ----------------------------------------------------------------------
# paper style (Computer Modern via matplotlib's bundled fonts)
# ----------------------------------------------------------------------
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
cnm  = [r"$\overline{u_1^2}/2k_t$", r"$\overline{u_2^2}/2k_t$", r"$\overline{u_3^2}/2k_t$"]

# ----------------------------------------------------------------------
# read one dump: returns (time, fractions[3], kt)
# ----------------------------------------------------------------------
def read_dump(fname):
    t = None
    pos=[]; posf=[]; u=[]; v=[]; w=[]
    with open(fname) as f:
        for line in f:
            if line.startswith("#"):
                m = re.search(r"time\s*=\s*([-\d.eE+]+)", line)
                if m: t = float(m.group(1))
                continue
            if not line.strip(): continue
            c = line.split()
            posf.append(float(c[1])); u.append(float(c[2]))
            v.append(float(c[3])); w.append(float(c[4]))
    posf = np.array(posf); u=np.array(u); v=np.array(v); w=np.array(w)
    # faces are the per-cell left faces (posf) plus the right domain boundary.
    # The homogeneousStrain domain is centred at 0, so the right boundary is
    # -posf[0] (= +L/2). Change this if your domain is not centred at 0.
    faces = np.append(posf, -posf[0])
    dx = np.diff(faces)
    L  = dx.sum()
    mu=(u*dx).sum()/L; mv=(v*dx).sum()/L; mw=(w*dx).sum()/L
    F=[u-mu, v-mv, w-mw]
    R=np.zeros((3,3))
    for a in range(3):
        for b in range(3):
            R[a,b]=(F[a]*F[b]*dx).sum()/L
    kt=0.5*np.trace(R)
    return t, np.array([R[i,i]/(2*kt) for i in range(3)]), kt

# ----------------------------------------------------------------------
# Level 0 benchmark for plane strain (exact RDT + LRR), normalised S=1
# ----------------------------------------------------------------------
I3=np.eye(3)
def planeA(): a=0.5; return np.diag([a,-a,0.0])
def production(A,R): return -(A@R+R@A.T)
def rapid_LRR(A,R):
    kt=0.5*np.trace(R); b=R/(2*kt)-I3/3; S=0.5*(A+A.T); W=0.5*(A-A.T)
    t2=b@S+S@b-(2/3)*np.trace(b@S)*I3; t3=W@b-b@W
    return 0.8*kt*S+1.75*kt*t2+1.31*kt*t3
def model_LRR(A,emax,n=200):
    te=np.linspace(0,emax,n)
    sol=solve_ivp(lambda t,y:(production(A,y.reshape(3,3))+rapid_LRR(A,y.reshape(3,3))).reshape(-1),
                  [0,emax],((2/3)*I3).reshape(-1),t_eval=te,rtol=1e-10,atol=1e-12)
    R=sol.y.reshape(3,3,-1); kt=0.5*(R[0,0]+R[1,1]+R[2,2])
    return te, np.array([R[i,i]/(2*kt) for i in range(3)])
def exact_rdt(A,emax,n=200,nmu=24,nphi=48):
    mu,wmu=leggauss(nmu); phi=(np.arange(nphi)+0.5)*2*np.pi/nphi
    dirs=[];w=[]
    for m,wm in zip(mu,wmu):
        st=np.sqrt(max(0,1-m*m))
        for p in phi: dirs.append([st*np.cos(p),st*np.sin(p),m]); w.append(0.5*wm/nphi)
    dirs=np.array(dirs);w=np.array(w);w/=w.sum();N=len(dirs)
    PHI0=I3[None]-np.einsum("ni,nj->nij",dirs,dirs)
    y0=np.concatenate([dirs.reshape(-1),PHI0.reshape(-1)])
    def rhs(t,y):
        K=y[:3*N].reshape(N,3);PH=y[3*N:].reshape(N,3,3)
        AtK=K@A;k2=np.einsum("ni,ni->n",K,K)
        M=-A[None]+2*np.einsum("ni,nj->nij",K,AtK)/k2[:,None,None]
        return np.concatenate([(-(K@A)).reshape(-1),
            (np.einsum("nik,nkj->nij",M,PH)+np.einsum("nik,njk->nij",PH,M)).reshape(-1)])
    te=np.linspace(0,emax,n)
    sol=solve_ivp(rhs,[0,emax],y0,t_eval=te,rtol=1e-9,atol=1e-11)
    PH=sol.y[3*N:].reshape(N,3,3,-1); R=np.einsum("n,nijt->ijt",w,PH)
    kt=0.5*(R[0,0]+R[1,1]+R[2,2])
    return te, np.array([R[i,i]/(2*kt) for i in range(3)])

# ----------------------------------------------------------------------
# gather dumps
# ----------------------------------------------------------------------
files = sorted(glob.glob(os.path.join(DUMP_DIR, "dmp_*.dat")))
for extra in ("odt_init.dat", "odt_end.dat"):
    p=os.path.join(DUMP_DIR, extra)
    if os.path.exists(p): files.append(p)

rows=[]
for f in files:
    try:
        t,fr,kt=read_dump(f)
        if t is not None: rows.append((t,fr,kt,os.path.basename(f)))
    except Exception as ex:
        print(f"skip {f}: {ex}")
rows.sort(key=lambda r:r[0])
e_odt = np.array([r[0]*SMAG for r in rows])
fr_odt= np.array([r[1] for r in rows])     # shape (Ndump,3)
kt_odt= np.array([r[2] for r in rows])
print(f"read {len(rows)} dumps, e in [{e_odt.min():.2f}, {e_odt.max():.2f}]")
emax=max(4.0, e_odt.max() if len(e_odt) else 4.0)

# ----------------------------------------------------------------------
# benchmark curves
# ----------------------------------------------------------------------
A=planeA()
te,fe=exact_rdt(A,emax)
tm,fm=model_LRR(A,emax)

# ----------------------------------------------------------------------
# plot
# ----------------------------------------------------------------------
fig,ax=plt.subplots(figsize=(5.4,3.9))
for i in range(3):
    ax.plot(te,fe[i],color=ccol[i],ls="-",lw=1.6)                 # exact RDT
    ax.plot(tm,fm[i],color=ccol[i],ls=(0,(6,2)),lw=1.3)           # LRR model
    ax.plot(e_odt,fr_odt[:,i],ls="none",marker="o",ms=3.4,        # ODT data
            mfc="none",mec=ccol[i],mew=0.9)
ax.set_xlabel(r"total strain $e = S\,t$")
ax.set_ylabel(r"component energy fraction")
ax.set_title(r"Strain-coupled ODT vs rapid-distortion theory (plane strain)")
ax.set_xlim(0,emax); ax.grid(alpha=0.25,lw=0.5)
leg_c=ax.legend([Line2D([0],[0],color=ccol[i],lw=2) for i in range(3)],cnm,
                title=r"component",loc="upper left",bbox_to_anchor=(1.02,1.0),
                frameon=False,handlelength=1.8)
ax.legend([Line2D([0],[0],color="k",ls="-",lw=1.6),
           Line2D([0],[0],color="k",ls=(0,(6,2)),lw=1.3),
           Line2D([0],[0],color="k",ls="none",marker="o",ms=4,mfc="none",mec="k")],
          [r"exact RDT",r"LRR closure",r"ODT (this work)"],
          title=r"source",loc="upper left",bbox_to_anchor=(1.02,0.55),
          frameon=False,handlelength=2.4)
ax.add_artist(leg_c)
for ext in ("pdf","png"):
    fig.savefig(f"{OUT}.{ext}",bbox_inches="tight")
print(f"wrote {OUT}.pdf / .png")

# also a TKE-amplification panel (diagnostic; not for the paper unless wanted)
fig2,ax2=plt.subplots(figsize=(5.0,3.5))
ax2.plot(e_odt,kt_odt/kt_odt[0] if len(kt_odt) else [],"o-",color="k",ms=3.4,lw=1.0)
ax2.set_xlabel(r"total strain $e=S\,t$"); ax2.set_ylabel(r"$k_t(e)/k_t(0)$")
ax2.set_title(r"Turbulent kinetic energy amplification"); ax2.grid(alpha=0.25,lw=0.5)
for ext in ("pdf","png"): fig2.savefig(f"{OUT}_kt.{ext}",bbox_inches="tight")
print(f"wrote {OUT}_kt.pdf / .png")
