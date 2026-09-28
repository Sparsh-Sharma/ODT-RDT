We are continuing the **ODT-RDT** project: strain-coupled One-Dimensional
Turbulence feeding an anisotropic upwash spectrum into Amiet leading-edge-noise
theory. Priority is the **JFM-2026-1781 revision** (Sharma & Kerstein); the
eddy-on-mean rebuild is a foundational side result.

Before doing anything, read, in this order:
1. `CLAUDE.md` (repo root) — entry point.
2. `docs/PROJECT_HANDBOOK.md` — the full architecture, branch map, HPC workflow,
   results, and decisions. Do not re-derive any of it.
3. The auto-memory index (loaded each session): read `odt-rdt-project-state`,
   `hpc-cara-caro-usage`, `eddymean-rebuild-verdict`, and the JFM-1781 set in
   full before running/building/campaigning ODT.
4. `docs/REVISION_MAP_jfm1781.md` — on `master` alongside the manuscript.

Branch map you must not get wrong (verify with `git worktree list` /
`git branch -avv`):
- **`master` = `origin/master` (a4d7dc6) = THE PAPER**, and the `ODT-RDT`
  working dir is on it: `manuscript/main_v4_overleaf.tex` (file of record),
  figures, `REVISION_MAP`, the `post/closure_bound/` pipeline, 128³ DNS. Do
  paper work here.
- **eddymean-rebuild (62dafbd)** = the eddy-on-mean rebuild solver + the T2
  verdict (`post/eddymean_results/`). All eddy-on-mean work stays on THIS
  branch — never commit it onto `master`.
- **`..\ODT-RDT-pm`** = a worktree detached at origin/master that also holds the
  six `fw_*_ensemble.npz` for any deep-strain follow-up.
- **`..\ODT-RDT-len` / `eddy-mean-forcing`** = a DIFFERENT active chat — do not
  touch it.
- **LEN_Extension** = Test-3 / clock machinery + `run/caro/` launchers +
  `docs/PLAN_paper2_test3.md`; the caro clone tracks it.
- **`master-eddymean-backup` (614a435)** = safety backup of a 2026-09-28 master
  divergence that was resolved (eddy-on-mean work had been committed onto master
  by mistake; master was reset to the paper). Deletable once confirmed unneeded.

HPC: all sims run on DLR **caro** (never locally). Connect via
`wsl bash -c "ssh -o BatchMode=yes caro '...'"` (login carologin6, VPN off-site,
SLURM account 2002095, partition `medium` = whole 256-core nodes). Pipe scripts
over stdin (`tr -d '\r' < file | ssh caro 'bash -s'`) — inline `$VAR`/loops
arrive empty. Rebuild + verify the binary after every pull. Launch 1024-rlz
cases with `run/caro/slrm_test3_node.sh`; give kvisc=1e-5 `fw_` cases
`--time=01:00:00`. Count completions via `odt_end.dat`, not `ls`. caro python:
`env PYTHONPATH= PYTHONNOUSERSITE=1 ~/anaconda3/envs/odt/bin/python3`.

Current state (2026-09-28):
- Paper: one-paper/three-part `main_v4_overleaf.tex` on `master` + ESM
  `supplementary.tex`; Alan Kerstein coauthor; awaiting his next review; some
  O-ledger items still open (see REVISION_MAP).
- Eddy-on-mean rebuild: T1 neutral, **T2 spectral NEGATIVE**, T3 cannot replace
  the relaxation machinery. Recommended decision (pending Alan): fold in only as
  a short foundational note (transition-from-laminar answers Alan's Rogallo
  objection), no spectral claim. Note to Alan drafted at
  `notes/email_alan_eddymean.html` (untracked).

Immediate options: (a) send Alan the eddy-on-mean note and get his read;
(b) on his verdict, draft the foundational-note paragraph into the manuscript;
(c) continue the open REVISION_MAP items; (d) the eddy-on-mean deep-strain
(e=2–3) follow-up, if wanted, using the ensembles in `..\ODT-RDT-pm`.

Conventions: commits NEVER carry a `Co-Authored-By:`/Claude line (author =
Sparsh Sharma only); a post-commit hook auto-pushes; figures follow
`post/closure_bound/figstyle_jfm.py`; writing is plain and measured, no
self-congratulation. Assess after each step before the next action.

Tell me where you want to start.
