# ODT & LEN — Claude Code migration package

This folder is a self-contained starter kit for moving the **ODT & LEN** research
project out of a claude.ai Project and into **Claude Code (desktop client)**.

It was assembled on 2026-08-23. It contains everything that could be exported from
the claude.ai Project (the reference literature as annotated notes, plus a full
context handover) together with the setup steps for wiring Claude Code to a local
project folder and to the two DLR HPC systems (**CARO** and **CARA**).

> **What is "LEN"?** The project is named *ODT & LEN*. "ODT" is unambiguously the
> One-Dimensional Turbulence model. "LEN" is your own shorthand — most likely the
> **Linear Eddy Model (LEM)**, which is ODT's direct ancestor (Kerstein 1991) and is
> referenced throughout the literature here. Wherever this kit says *LEN*, read it as
> whatever code/model that label means to you, and correct the one line in `CLAUDE.md`
> if the guess is wrong.

---

## Read in this order

1. **`HANDOVER.md`** — the full project context: what the work is, the physics, the
   code lineage, the papers, and the open threads. Read this first, and paste its
   summary into your first Claude Code session.
2. **`docs/claude-code-migration-guide.md`** — step-by-step: install Claude Code,
   choose and connect the project folder, and connect CARO + CARA over SSH.
3. **`docs/recommended-structure.md`** — the proposed repository layout for the code,
   cases, data, and notes.
4. **`docs/hpc-cheatsheet.md`** — CARO and CARA quick reference (SLURM, modules,
   file transfer).
5. **`papers/`** — annotated bibliography (`PAPERS.md`) and a note on bringing the
   original PDFs across (`README.md`).
6. **`prompts/first-session-prompt.md`** — a ready-to-paste kickoff prompt for the
   first Claude Code session.

## The three things only you can supply

This kit is complete except for material that never lived in the claude.ai Project and
therefore could not be exported from it:

- **The 6 original PDFs** — they live in the claude.ai Project. See `papers/README.md`.
- **The ODT/LEN source code** — lives on your machine and/or on CARO/CARA.
- **The simulation data and your research notes** — on the clusters / your disk.

`docs/recommended-structure.md` shows exactly where each of these lands in the new
project.
