#!/usr/bin/env python3
r"""
Summarise the resolution sweep produced by run_sweep.sh.

For each <root>/hs_*/data, reads the first (e=0) and last dump, computes the
TKE amplification kt(e_last)/kt(0), and compares it to the LRR moment model
integrated to the SAME e_last (so the deficit is matched in strain, not against
a fixed e=4 value). Also reports the final component-energy fractions.

numpy only (no scipy). Usage:  python3 analyze_sweep.py [root]   (default ../data)
"""
import os, glob, re, sys
import numpy as np

ROOT = sys.argv[1] if len(sys.argv) > 1 else "../data"
I3 = np.eye(3)
A  = np.diag([0.5, -0.5, 0.0])          # plane strain, S = 1  (e = t)

# ---- LRR moment model integrated with RK4 (numpy) to arbitrary e ----
def lrr_rhs(R):
    P = -(A @ R + R @ A.T)
    kt = 0.5*np.trace(R); b = R/(2*kt) - I3/3.0
    S = 0.5*(A + A.T)                    # W = 0 (irrotational)
    t2 = b@S + S@b - (2.0/3.0)*np.trace(b@S)*I3
    return P + 0.8*kt*S + 1.75*kt*t2
def lrr_to(e, de=0.005):
    R = (2.0/3.0)*I3; n = max(1, int(round(e/de))); h = e/n
    for _ in range(n):
        k1=lrr_rhs(R); k2=lrr_rhs(R+0.5*h*k1); k3=lrr_rhs(R+0.5*h*k2); k4=lrr_rhs(R+h*k3)
        R = R + (h/6.0)*(k1+2*k2+2*k3+k4)
    kt=0.5*np.trace(R)
    return kt/1.0, np.array([R[i,i]/(2*kt) for i in range(3)])   # kt(0)=1

# ---- read one ODT dump ----
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
    posf=np.array(posf); u=np.array(u); v=np.array(v); w=np.array(w)
    faces=np.append(posf, -posf[0]); dx=np.diff(faces); L=dx.sum()
    mu=(u*dx).sum()/L; mv=(v*dx).sum()/L; mw=(w*dx).sum()/L
    F=[u-mu, v-mv, w-mw]; R=np.zeros((3,3))
    for a in range(3):
        for b in range(3): R[a,b]=(F[a]*F[b]*dx).sum()/L
    kt=0.5*np.trace(R)
    return t, kt, np.array([R[i,i]/(2*kt) for i in range(3)]), len(dx)

def dumps_in(d):
    fs=glob.glob(os.path.join(d,"data_*","dmp_*.dat"))
    return fs if fs else glob.glob(os.path.join(d,"dmp_*.dat"))

def has_dumps(case_dir):
    return bool(dumps_in(os.path.join(case_dir,"data")))

# Accept either a single case directory (ROOT/data/... holds dumps) or a root
# whose immediate children are case directories (hs_*, homogeneousStrain*, ...).
if has_dumps(ROOT):
    cases=[ROOT]
else:
    cases=sorted(d for d in glob.glob(os.path.join(ROOT,"*"))
                 if os.path.isdir(d) and has_dumps(d))
if not cases:
    print(f"no case directories with dumps found at/under {ROOT}")
    print("point me at a case dir (…/homogeneousStrain1b) or a root (…/data)")
    sys.exit(1)

print(f"{'config':12s}{'e_last':>7s}{'ngrd':>7s}{'kt0':>8s}{'ODT amp':>9s}"
      f"{'LRR amp':>9s}{'deficit':>9s}   {'f11':>7s}{'f22':>7s}{'f33':>7s}")
print("-"*92)
for c in cases:
    label=os.path.basename(c).replace("hs_","")
    fs=dumps_in(os.path.join(c,"data"))
    if not fs: print(f"{label:12s}  (no dumps)"); continue
    rows=[]
    for f in fs:
        try:
            t,kt,fr,ng=read_dump(f)
            if t is not None: rows.append((t,kt,fr,ng))
        except Exception: pass
    rows.sort(key=lambda r:r[0])
    _,kt0,_,_ = rows[0]; eN,ktN,frN,ngN = rows[-1]
    amp = ktN/kt0
    lrr_amp,_ = lrr_to(eN)
    deficit = 100.0*(amp-lrr_amp)/lrr_amp
    print(f"{label:12s}{eN:7.2f}{ngN:7d}{kt0:8.3f}{amp:9.3f}{lrr_amp:9.3f}"
          f"{deficit:8.1f}%   {frN[0]:7.4f}{frN[1]:7.4f}{frN[2]:7.4f}")
print("-"*92)
print("deficit = (ODT amp - LRR amp)/LRR amp, both at e_last (matched strain).")
print("Flat deficit across configs => NOT a spatial-resolution effect.")