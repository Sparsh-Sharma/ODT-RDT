# CLAUDE.md — ODT & LEN project

> This file is read automatically by Claude Code at the start of every session in this
> repository. Keep it short, current, and factual. Put long-form context in
> `HANDOVER.md` and point to it. Move this file to the **root of your actual project
> repository** once you have set it up (see `docs/recommended-structure.md`).

## What this project is

Reduced-order **stochastic turbulence modelling** for **aeroacoustics** — using the
**One-Dimensional Turbulence (ODT)** model (and the related **LEN**/Linear-Eddy-Model
lineage) as a full-scale-resolving, 1-D surrogate for DNS to generate turbulent
velocity fields, then feeding those fields into an acoustic analogy to predict
**far-field sound radiation** (jet noise, low-Mach-number jets to start).

Owner: Dr. Sparsh Sharma (DLR Braunschweig). Long-running collaboration with the
BTU Cottbus-Senftenberg ODT group (Chair of Numerical Fluid and Gas Dynamics).

> ⚠️ Confirm the meaning of **LEN** and edit this line if needed. Assumed here to be
> the Linear Eddy Model (LEM), Kerstein 1991 — ODT's predecessor.

## The model in one paragraph

ODT evolves velocity/scalar profiles on a 1-D line of sight through a 3-D turbulent
flow. Molecular diffusion is solved on the line; turbulent advection is replaced by a
stochastic sequence of instantaneous **eddy events**, each a **triplet map** (compress
an interval `[x0, x0+l]` to a third, place three copies, invert the middle) plus a
kernel operation (`J`, `K`) that redistributes energy between velocity components to
model pressure scrambling / return-to-isotropy. Eddies are sampled from a rate
distribution set by available kinetic energy. Three model parameters: **C** (eddy-rate
/ turbulence intensity), **Z** (viscous-penalty cutoff for small eddies), and a
large-eddy suppression parameter (**βLES** / elapsed-time method). Modern
implementations use an adaptive Lagrangian finite-volume mesh.

## Reference code

- **BYUignite/ODT** — David Lignell's open-source, object-oriented **C++** ODT code
  (Stephens & Lignell, SoftwareX 2021). MIT licensed. Build: **CMake ≥ 3.12**,
  **Cantera**, **yaml-cpp**, optional Doxygen. Input via `input.yaml`. Runs many
  independent **realizations** (embarrassingly parallel) via MPI/SLURM
  (`slrmJob.sh`, `slrmJob_array.sh`). Post-processing in Python.
- A BTU-maintained variant/branch of the ODT code is the likely working codebase.
  **Record its git remote and build quirks here once confirmed.**

## HPC

Runs on the two DLR clusters, both AMD EPYC, both SLURM, both usable by all DLR
institutes:

- **CARO** — Göttingen (GWDG / Univ. Göttingen data centre). ~1,364 nodes, 64 cores/node.
- **CARA** — TU Dresden (ZIH / Lehmann centre). ~2,280 nodes (mostly 32 cores/node) + A100 GPU nodes.

Login hostnames, project/account IDs, and module names are institute-specific — see
`docs/hpc-cheatsheet.md` and fill in the placeholders.

## Conventions & working agreements

- ODT produces **ensembles of realizations**; nothing is meaningful until averaged.
  Always keep the number of realizations and RNG seeds with any result.
- Keep `C`, `Z`, `βLES`, Reynolds/Mach numbers, and realization count attached to every
  case and every figure.
- Do **not** run heavy simulations inside Claude Code's local sandbox — Claude Code
  edits code and drives job submission; the compute happens on CARO/CARA via SLURM.
- Big data (`data/`, `runtime/`, raw dumps) stays on the cluster / out of git — see
  `.gitignore`.

## Owner preferences (from prior sessions)

- Wants deep technical rigor and physics-first explanations; push back on imprecise
  framings rather than smoothing them over.
- Prefers native/first-principles derivations over hand-waving.

## Pointers

- Full context: `HANDOVER.md`
- Setup & HPC wiring: `docs/claude-code-migration-guide.md`, `docs/hpc-cheatsheet.md`
- Literature: `papers/PAPERS.md`
