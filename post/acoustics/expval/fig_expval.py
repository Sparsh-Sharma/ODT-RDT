#!/usr/bin/env python3
"""Four-way validation figure: measured baseline LE-noise spectra vs
(1) classical Amiet+vK, (2) Amiet+Liepmann, (3) Amiet+standard ODT,
(4) Amiet+fixed ODT (size-capped relaxation clock), per case.

Usage: fig_expval.py [case ...]   (default: all cases with a csv file)
Writes fig_expval_<case>.{pdf,png} + metrics_expval.txt (mean and rms
error over the stated comparison band, per model).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..",
                                                 "closure_bound")))
import chain  # noqa: E402
from cases import CASES, load_csv  # noqa: E402
from figstyle_jfm import COL, FULL, GOLD, RED, plt, save  # noqa: E402

MODELS = [("vk", "Amiet + von Kármán", "k", (0, (5, 2))),
          ("liepmann", "Amiet + Liepmann", "0.45", (0, (1, 1.2))),
          ("odt_std", "Amiet + ODT (standard)", COL["odt"], "-"),
          ("odt_fix", "Amiet + ODT (fixed)", RED, "-")]


def have_fix(tag):
    return os.path.exists(os.path.join(
        chain.LENOISE, f"rdt_family_fit_table_{tag}_RCS1.txt"))


def run_case(key, ax=None):
    c = CASES[key]
    fexp, spl_exp = load_csv(c["csv"])
    f = np.geomspace(max(fexp.min(), 80.0), fexp.max(), 48)
    out, metrics = {}, {}
    for m, lab, col, ls in MODELS:
        if m == "odt_fix" and not have_fix(c["fit_tag"]):
            continue
        kw = dict(t_over_c=c["t_over_c"], fit_tag=c["fit_tag"])
        if m.startswith("odt"):
            kw["e_eff"] = c["e_ctr"]
        spl = chain.predict_spl(f, c["U"], c["w2"], c["Lam"], c["chord"],
                                c["span"], c["obs"], m, **kw) \
            + c["slc_db"]
        out[m] = spl
        # metrics on the stated band
        b = c["band"]
        sel = (fexp >= b[0]) & (fexp <= b[1])
        err = np.interp(fexp[sel], f, spl) - spl_exp[sel]
        metrics[m] = (float(err.mean()),
                      float(np.sqrt((err ** 2).mean())),
                      float(err.std()))
    if ax is not None:
        b = c["band"]
        inb = (fexp >= b[0]) & (fexp <= b[1])
        ax.plot(fexp[inb], spl_exp[inb], "o", ms=3, mfc="none", mec="k",
                mew=0.7, label="measured")
        ax.plot(fexp[~inb], spl_exp[~inb], "o", ms=3, mfc="none",
                mec="0.65", mew=0.6)
        for m, lab, col, ls in MODELS:
            if m in out:
                ax.plot(f, out[m], color=col, ls=ls, lw=1.1, label=lab)
        # e_eff sensitivity band for the ODT curves
        if "odt_std" in out and c["e_hi"] > 0:
            lo = chain.predict_spl(f, c["U"], c["w2"], c["Lam"],
                                   c["chord"], c["span"], c["obs"],
                                   "odt_std", t_over_c=c["t_over_c"],
                                   fit_tag=c["fit_tag"],
                                   e_eff=c["e_lo"]) + c["slc_db"]
            hi = chain.predict_spl(f, c["U"], c["w2"], c["Lam"],
                                   c["chord"], c["span"], c["obs"],
                                   "odt_std", t_over_c=c["t_over_c"],
                                   fit_tag=c["fit_tag"],
                                   e_eff=c["e_hi"]) + c["slc_db"]
            ax.fill_between(f, lo, hi, color=COL["odt"], alpha=0.15,
                            lw=0)
        ax.set_xscale("log")
        ax.set_xlabel(r"$f$ [Hz]")
        ax.set_ylabel(r"SPL [dB Hz$^{-1}$ re $20\,\mu$Pa]")
        ax.set_title(c["name"], fontsize=8)
    return metrics


def main():
    keys = sys.argv[1:] or [k for k in CASES
                            if os.path.exists(CASES[k]["csv"])]
    lines = []
    for k in keys:
        fig, ax = plt.subplots(figsize=(0.62 * FULL, 2.6))
        met = run_case(k, ax)
        ax.legend(fontsize=6, loc="lower left")
        fig.tight_layout(pad=0.4)
        save(fig, os.path.join(HERE, f"fig_expval_{k}"))
        plt.close(fig)
        lines.append(f"== {CASES[k]['name']}  band {CASES[k]['band']}")
        for m, (bias, rms, shape) in met.items():
            lines.append(f"   {m:9s} bias {bias:+6.2f} dB   rms "
                         f"{rms:5.2f} dB   shape {shape:5.2f} dB")
    txt = "\n".join(lines)
    open(os.path.join(HERE, "metrics_expval.txt"), "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
