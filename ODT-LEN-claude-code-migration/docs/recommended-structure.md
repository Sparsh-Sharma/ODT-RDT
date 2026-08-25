# Recommended repository structure

A layout that keeps code, cases, notes, and HPC glue in one Claude Code workspace while
keeping large simulation data **out** of git (it lives on CARO/CARA).

```
odt-len/
├── CLAUDE.md                  # project memory (Claude Code reads at root) — from this kit
├── HANDOVER.md                # full context — from this kit
├── README.md                  # optional: your own top-level readme
├── .gitignore                 # from this kit (excludes build/, data/, PDFs, secrets)
│
├── code/
│   ├── odt/                   # BYUignite/ODT clone or your BTU/DLR fork (its own git)
│   └── len/                   # LEN/LEM code, if a separate codebase  ← confirm what LEN is
│
├── cases/                     # case DEFINITIONS only (small, versioned)
│   └── <case-name>/
│       ├── input.yaml         # ODT parameters (C, Z, βLES, Re, Ma, nRlz, seeds…)
│       └── notes.md           # what this case is, expected outcome, status
│
├── acoustics/                 # the SPL / far-field pipeline (Medina et al. 2023)
│   ├── src/                   # velocity-field → spectral sources → SPL
│   └── validation/            # e.g. Viswanathan (2007) Ma=0.4 jet comparison
│
├── post/                      # post-processing scripts (Python) — versioned
│
├── results/                   # SMALL post-processed outputs + figures pulled back
│   └── <case-name>/           # (raw dumps stay on the cluster, not here)
│
├── notes/                     # research notes, derivations, drafts, paper-in-progress
│
├── hpc/
│   ├── ssh_config.sample      # caro/cara Host blocks (no secrets)
│   ├── slurm/                 # job scripts (adapt from ODT's slrmJob*.sh)
│   └── env/                   # module-load snippets per cluster
│
├── papers/                    # literature — from this kit (+ the 6 PDFs)
│   ├── PAPERS.md
│   └── *.pdf
│
└── docs/                      # this kit's docs
```

## What lives where — the rule

| Thing | In git / local folder? | Where it really lives |
|---|---|---|
| ODT/LEN source code | ✅ (as submodule or nested clone) | its own git remote |
| Case definitions (`input.yaml`) | ✅ small, versioned | here |
| Post-processing + acoustics scripts | ✅ versioned | here |
| Research notes / drafts | ✅ versioned | here |
| **Raw simulation data, `runtime/`, dumps** | ❌ `.gitignore` | **CARO/CARA scratch** |
| Post-processed results + figures | ✅ only the small ones | `results/` + cluster |
| Cluster hostnames / account IDs / keys | ❌ never in git | `~/.ssh/config`, local env |

## Notes
- Treat `code/odt/` as an upstream you rebase onto; keep **your** changes as a fork/branch
  and record the remote + branch in `CLAUDE.md`.
- If LEN is a separate model/codebase, mirror this structure under `code/len/` and give it
  its own short `CLAUDE.md` note. If "LEN" just refers to the LEM lineage inside ODT,
  delete `code/len/` and say so in `CLAUDE.md`.
- Keep every result reproducible: the `input.yaml`, the code commit hash, the realization
  count, and the seeds together. A one-line `provenance.txt` per results folder is enough.
