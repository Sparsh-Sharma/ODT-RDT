#!/usr/bin/env python3
r"""
Summarise the kvisc0 sweep (run_kvisc_sweep.sh) for the eddies-ON + dilatation-ON
runs. For each ../data/hs2_kv*/ it reads the first and last dump, computes the
spectral centroid migration kbar(e_last)/kbar(0), and compares it to the linear-
RDT compression exp(-A22 e_last). This locates the Reynolds-number window:

   centroid << RDT  -> too viscous (dissipation kills the compressed energy)
   centroid ~  RDT  -> eddies negligible (behaves ~linear)
   centroid >  RDT  -> cascade adds transfer beyond compression  [the target]

It also reports ngrd_end (resolution) and the dumps count (did it complete).

Usage:  python3 analyze_kvisc_sweep.py [root]      (root defaults to ../data)
"""
import os, glob, re, sys
import numpy as np

ROOT     = sys.argv[1] if len(sys.argv) > 1 else "../data"
SMAG     = 1.0
A22      = -0.5
NUNIFORM = 2048
CENT_THRESH = 1e-3

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

def centroid(posf, fields, Nu):
    x0=posf[0]; L=-2.0*x0
    faces=np.append(posf,-x0); xc=0.5*(faces[:-1]+faces[1:])
    xu=x0+(np.arange(Nu)+0.5)*(L/Nu)
    k=2.0*np.pi*np.fft.rfftfreq(Nu,d=L/Nu)
    cs=[]
    for f in fields:
        fu=np.interp(xu,xc,f,period=L); fu=fu-fu.mean()
        E=2.0*np.abs(np.fft.rfft(fu)/Nu)**2
        E=E[1:]; kk=k[1:]
        m=E>CENT_THRESH*E.max()
        cs.append((kk[m]*E[m]).sum()/E[m].sum() if m.any() else np.nan)
    return cs

def dumps_in(d):
    fs=glob.glob(os.path.join(d,"data_*","dmp_*.dat"))
    return fs if fs else glob.glob(os.path.join(d,"dmp_*.dat"))

cases=sorted(glob.glob(os.path.join(ROOT,"hs2_*")))
if not cases:
    print(f"no hs2_* directories under {ROOT}"); sys.exit(1)

print(f"{'case':10s}{'e_last':>7s}{'ngrd':>7s}{'dumps':>7s}"
      f"{'cbar1/0':>9s}{'cbar2/0':>9s}{'cbar3/0':>9s}{'RDT':>8s}{'verdict':>14s}")
print("-"*86)
for c in cases:
    label=os.path.basename(c).replace("hs2_","")
    fs=dumps_in(os.path.join(c,"data"))
    nd=len(fs)
    if nd==0:
        print(f"{label:10s}  (no dumps -- crashed/early-exit)"); continue
    rows=[]
    for f in fs:
        t,posf,u,v,w=read_dump(f)
        if t is not None: rows.append((t,posf,u,v,w))
    rows.sort(key=lambda r:r[0])
    t0,p0,u0,v0,w0=rows[0]; tN,pN,uN,vN,wN=rows[-1]
    eN=tN*SMAG
    c0=centroid(p0,[u0,v0,w0],NUNIFORM)
    cN=centroid(pN,[uN,vN,wN],NUNIFORM)
    rat=[cN[i]/c0[i] for i in range(3)]
    rdt=np.exp(-A22*eN)
    avg=np.nanmean(rat)
    if   avg < 0.6*rdt: verdict="too viscous"
    elif avg > 1.05*rdt: verdict="cascade>RDT *"
    else:               verdict="~ RDT"
    print(f"{label:10s}{eN:7.2f}{len(pN):7d}{nd:7d}"
          f"{rat[0]:9.3f}{rat[1]:9.3f}{rat[2]:9.3f}{rdt:8.2f}{verdict:>14s}")
print("-"*86)
print("cbar_i/0 = spectral-centroid migration; RDT = exp(-A22 e_last) (pure compression).")
print("Target window: cbar climbing toward/above RDT with the run completing (~40 dumps)")
print("and ngrd not collapsing. Then inspect that case's spectrum with spectrum_diagnostic.py.")
