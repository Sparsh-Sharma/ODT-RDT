#!/usr/bin/env python3
"""Validation-case registry: published LE-noise experiments with fully
documented configuration and digitised baseline spectra (see
../expdata/CASES.md and the extraction scripts there).

Each entry: one velocity of one experiment.
  w2      : mean-square upwash u'_w^2 [(m/s)^2] (vertical component
            where reported, else Tu*U)
  Lam     : streamwise (longitudinal) integral length scale [m]
  obs     : observer (x1, x2, x3) [m] rel. LE/mid-span, x3 = overhead
  t_over_c: thickness ratio for the Gershfeld response admittance
            (None = flat plate)
  fit_tag : ODT ensemble family for the representation (rapidity)
  e_ctr, e_lo, e_hi : effective free-distortion strain (potential-flow
            stagnation-line accumulation truncated at the eddy scale;
            centre x*=Lambda/2, band x* in [Lambda, Lambda/4])
  csv     : digitised measured baseline spectrum (f_Hz, SPL_dB)
  band    : [fmin, fmax] Hz for quantitative metrics (above facility
            contamination, below background flags and the thickness-
            breakdown frequency U/t_A)
  slc_db  : shear-layer amplitude correction added to predictions
            (theory -> as-measured convention of the experiment)
"""
import os

import numpy as np

from chain import e_eff_ellipse

EXP = os.path.normpath(os.path.join(os.path.dirname(
    os.path.abspath(__file__)), "..", "expdata"))


def _e3(chord, thick, lam):
    return (e_eff_ellipse(chord, thick, lam / 2),
            e_eff_ellipse(chord, thick, lam),
            e_eff_ellipse(chord, thick, lam / 4))


def paterson_amiet(U, tu_w, lam, band, slc):
    e_ctr, e_lo, e_hi = _e3(0.23, 0.0276, lam)
    return dict(
        name=f"Paterson-Amiet 1976, NACA0012, U={U:.0f} m/s",
        U=U, w2=(tu_w * U) ** 2, Lam=lam, chord=0.23, span=0.53,
        obs=(0.0, 0.0, 2.25), t_over_c=0.12, fit_tag="S20",
        e_ctr=e_ctr, e_lo=e_lo, e_hi=e_hi,
        csv=os.path.join(EXP, f"PatersonAmiet_fig13_U{U:.0f}.csv"),
        band=band, slc_db=slc)


def bampanis22(U):
    e_ctr, e_lo, e_hi = _e3(0.10, 0.003, 0.009)
    return dict(
        name=f"Bampanis 2022 (ECL), flat plate, U={U:.0f} m/s",
        U=U, w2=(0.045 * U) ** 2, Lam=0.009, chord=0.10, span=0.30,
        obs=(0.0, 0.0, 1.25), t_over_c=None, fit_tag="S20",
        e_ctr=e_ctr, e_lo=e_lo, e_hi=e_hi,
        csv=os.path.join(EXP, f"Bampanis2022_fig3b_U{U:.0f}.csv"),
        band=[500.0, 8000.0], slc_db=0.0)


def bampanis19(U):
    c = bampanis22(U)
    c["name"] = f"Bampanis 2019 (ECL), flat plate, U={U:.0f} m/s"
    c["csv"] = os.path.join(EXP, f"Bampanis2019_fig3a_U{U:.0f}.csv")
    return c


def narayanan15():
    e_ctr, e_lo, e_hi = _e3(0.15, 0.002, 0.006)
    return dict(
        name="Narayanan 2015 (ISVR), flat plate, U=60 m/s",
        U=60.0, w2=(0.025 * 60) ** 2, Lam=0.006, chord=0.15, span=0.45,
        obs=(0.0, 0.0, 1.2), t_over_c=None, fit_tag="S20",
        e_ctr=e_ctr, e_lo=e_lo, e_hi=e_hi,
        csv=os.path.join(EXP, "Narayanan2015_fig4_baseline60.csv"),
        band=[500.0, 8000.0], slc_db=0.0)


CASES = {
    "PA40": paterson_amiet(40.0, 0.0453, 0.0302, [250, 1400], 0.15),
    "PA60": paterson_amiet(60.0, 0.0392, 0.0301, [250, 1500], 0.25),
    "PA90": paterson_amiet(90.0, 0.0482, 0.0294, [250, 2100], 0.35),
    "BA22_19": bampanis22(19.0),
    "BA22_27": bampanis22(27.0),
    "BA22_32": bampanis22(32.0),
    "BA19_32": bampanis19(32.0),
    "NA15_60": narayanan15(),
}


def load_csv(path):
    d = np.genfromtxt(path, delimiter=",", names=True)
    return d["f_Hz"], d["SPL_dB"]
