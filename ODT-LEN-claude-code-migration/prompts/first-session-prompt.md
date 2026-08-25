# First Claude Code session — kickoff prompt

Paste this into your first Claude Code session in the `odt-len/` folder (edit the
bracketed bits first). `CLAUDE.md` at the repo root already gives Claude the standing
context; this just sets the first task.

---

```
We're continuing my ODT & LEN research project, now in Claude Code. Read CLAUDE.md,
HANDOVER.md, and papers/PAPERS.md first for context — this is reduced-order stochastic
turbulence modelling (One-Dimensional Turbulence) coupled to an acoustic analogy for
far-field jet-noise prediction, in collaboration with the BTU Cottbus ODT group.

Environment:
- The ODT code lives at code/odt/ (remote: <git remote>, branch: <branch>).
- I run realizations on DLR CARO and CARA (both SLURM). SSH aliases `caro` and `cara`
  are set up. Big data stays on the clusters; this repo holds code, cases, scripts, notes.

Before writing anything, confirm you can:
1. Summarize, in your own words, the ODT eddy-event mechanism (triplet map + kernel)
   and the roles of the C, Z, and βLES parameters — so I know the context loaded.
2. Locate the SLURM run scripts in the code and the acoustics/SPL pipeline.

Then let's start on: <YOUR FIRST TASK — e.g. "reproduce the Ma=0.4, Re≈553k jet SPL
case from Medina Méndez et al. 2023 and set up a clean cases/ + hpc/slurm/ workflow for
it on CARA">.

Please push back on any imprecise framing rather than smoothing it over, and keep
explanations physics-first.
```

---

## Good first tasks to choose from
- Reproduce the **PAMM 2023** validation case (`Ma_j=0.4`, `Re≈553k`) end-to-end and
  lock down a clean `cases/` + `hpc/slurm/` workflow.
- Audit the acoustics pipeline against the paper's Eqs. (14)–(17) and write tests.
- Set up the round-jet velocity-statistics case from **Sharma, Klein & Schmidt (2022)**
  as the upstream input to the acoustics.
- Confirm what **LEN** is in this project and, if it's a separate codebase, get it
  building and documented under `code/len/`.
- Extend the framework toward **variable-density / heated jets** (Ashurst–Kerstein 2005
  machinery) — scope the changes needed in the code and the source terms.
