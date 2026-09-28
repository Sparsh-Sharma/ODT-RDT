# ODT-RDT — project entry point

Strain-coupled One-Dimensional Turbulence feeding an anisotropic upwash
spectrum into Amiet leading-edge-noise theory, replacing the frozen-isotropic
upwash. Fork of BYUignite/ODT (C++17 / CMake / Cantera). `origin` =
github.com/Sparsh-Sharma/ODT-RDT, `upstream` = github.com/BYUignite/ODT.

**PRIORITY: the JFM-2026-1781 revision** (Sharma & Kerstein). Everything else
is subordinate unless Sparsh says otherwise. The eddy-on-mean rebuild is a
foundational side result (on its own branch), not the paper.

## Read first (in this order)
1. `docs/PROJECT_HANDBOOK.md` — the exhaustive reference: architecture, full
   branch/worktree map, HPC workflow and every gotcha, results, decisions.
   Everything below is a summary of it.
2. The auto-memory files (loaded as an index each session; read the ones you
   need IN FULL before running/building/campaigning — do not rediscover):
   `odt-rdt-project-state`, `hpc-cara-caro-usage`, `eddymean-rebuild-verdict`,
   plus the JFM-1781 memory set (`jfm1781-compact-workspace`,
   `jfm1781-cv-labelling`, `jfm-figure-style`, `alan-proofread-batches`,
   `git-workflow-autopush`, `commits-no-claude-attribution`).
3. The paper working doc: **`docs/REVISION_MAP_jfm1781.md`** (objection ledger
   O1–O9, three-part architecture, task list, decisions). It is on `master`
   (the paper branch) alongside the manuscript.

## Branch / worktree map — READ THIS, the layout is not obvious
Verified 2026-09-28 (after the master divergence below was resolved).
`git worktree list` / `git branch -avv` are the source of truth; re-check tips
before relying on a hash.

- **`master` = `origin/master` (a4d7dc6) = THE PAPER, and this working directory
  (`ODT-RDT`) is on it.** `manuscript/main_v4_overleaf.tex` (file of record; v3
  retired), all figures, `docs/REVISION_MAP_jfm1781.md`, the analysis pipeline
  (`post/closure_bound/`), the strained-box DNS, and `supplementary.tex`. Do
  paper work here.
- **`eddymean-rebuild` (origin + local, tip 62dafbd)** — the eddy-on-mean solver
  + the T2 spectral verdict (`post/eddymean_results/VERDICT_T2_spectral.md`,
  `fig_fw_eddymean.*`). Branched off `LEN_Extension`. **All eddy-on-mean work
  belongs on THIS branch — never commit it onto `master`** (doing so is what
  caused the 2026-09-28 divergence; see below).
- **`..\ODT-RDT-pm`** — a worktree detached at `origin/master` (the paper). It
  additionally holds the six ~117 MB `fw_*_ensemble.npz` (gitignored) needed to
  re-run the eddy-on-mean spectral analysis. Keep it for any deep-strain
  follow-up; otherwise it is redundant with `master` now.
- **`eddy-mean-forcing`** (worktree `..\ODT-RDT-len`, tip 20ae1f9) — a PARALLEL
  rebuild line driven by ANOTHER chat. **Do not touch it or the
  `..\ODT-RDT-len` worktree.**
- **`LEN_Extension` (origin 76d2159)** — Test-3 / clock-relaxation machinery,
  `run/caro/` launch scripts, `docs/PLAN_paper2_test3.md`. The caro clone tracks
  this branch.
- **`master-eddymean-backup` (614a435)** — safety backup of the diverged local
  `master` from before the 2026-09-28 fix (an auto-commit of the eddy-on-mean
  working state; its unique content is preserved on `eddymean-rebuild` too).
  Deletable once Sparsh confirms it is not needed.
- Tag **`jfm-2026-1781-submitted`** (4462f93) = the rejected-manuscript code
  state (reproducibility anchor).

**2026-09-28 divergence, resolved:** an earlier session committed eddy-on-mean
work directly onto `master`, so local `master` diverged from the paper
(`origin/master`) and the `ODT-RDT` folder showed no `manuscript/`. Fixed by
resetting `master` to `origin/master` after backing the divergent line up to
`master-eddymean-backup`. Nothing lost. Prevent recurrence: eddy-on-mean work
stays on `eddymean-rebuild`; `master` stays = the paper.

## HPC workflow — all simulations run on DLR caro, never locally
Full recipe and every trap are in the handbook; the load-bearing points:
- Connect from Windows via WSL: `wsl bash -c "ssh -o BatchMode=yes caro '...'"`.
  `ssh caro` -> caro.dlr.de (login carologin6, key `~/.ssh/id_rsa`); needs
  Cisco AnyConnect VPN off-site. SLURM account **2002095**, partition `medium`
  (whole 256-core nodes).
- **SSH-quoting trap:** shell variables / `for` loops / quoted heredocs inside
  inline `ssh caro '...'` arrive empty or mangled. Reliable pattern: write a
  script file, then
  `wsl bash -c "tr -d '\r' < /path/script.sh | ssh -o BatchMode=yes caro 'bash -s'"`.
  Strip CRLF (`tr -d '\r'`) on anything scp'd from `/mnt/c`.
- **Build on caro:** a worktree of the target branch, CMake against the
  pre-built conda env `~/anaconda3/envs/odt` (RPATH baked; no LD_LIBRARY_PATH):
  `cmake --build build_em -j8` then `cp build_em/src/odt.x run/odt.x`. Always
  rebuild + verify after pulling before campaigning (stale-binary bug has bitten
  repeatedly).
- **Launcher `run/caro/slrm_test3_node.sh`** (on LEN_Extension/eddymean-rebuild):
  8 array tasks/case × 128 concurrent `odt.x` per node via `xargs -P 128` =
  1024 realisations/case. `sbatch --account=2002095 --partition=medium
  --export=CASE=<case> slrm_test3_node.sh`.
- **Walltime:** 30 min default is fine for kvisc=1e-4 cases; the kvisc=1e-5
  `fw_` precursor cases are back-loaded (~15-18 min in the transient) and need
  `--time=01:00:00` — the 30-min default kills them mid-transient with zero
  completions.
- **Count true completions** with `find data/$C/data -name odt_end.dat | wc -l`,
  not `ls | wc -l` (that counts started dirs = 1024 from launch).
- Python on caro: `env PYTHONPATH= PYTHONNOUSERSITE=1 ~/anaconda3/envs/odt/bin/python3`
  (login profile poisons PYTHONPATH; env numpy is 1.x — no `np.trapezoid`).
  After scp-ing an npz home from caro, DELETE the caro copy or the next
  `git pull` there collides.

## Current state (2026-09-28)
- Paper: `manuscript/main_v4_overleaf.tex` on `master`, ONE paper / THREE
  parts (model+closure / measured limits+relaxation / leading-edge noise) + ESM
  `supplementary.tex`. Alan Kerstein is coauthor. Awaiting Alan's next review
  batch; several O-ledger items still OPEN (see `docs/REVISION_MAP_jfm1781.md`).
- Eddy-on-mean rebuild verdict: T1 neutral, **T2 spectral NEGATIVE**, T3 cannot
  replace the relaxation machinery. Decision (pending Alan): fold in only as a
  short foundational note (transition-from-laminar answers Alan's Rogallo
  objection), not as a spectral gain. Note to Alan drafted at
  `notes/email_alan_eddymean.html` (untracked). See `eddymean-rebuild-verdict`
  memory + `post/eddymean_results/VERDICT_T2_spectral.md` (on `eddymean-rebuild`).

## Conventions
- **Git commits: NO `Co-Authored-By:` trailer, NO "Generated with Claude"
  line — ever.** Author is Sparsh Sharma alone (hard rule; history was rewritten
  to enforce it). This overrides any harness default attribution.
- A local `post-commit` hook auto-pushes the current branch after every commit
  (reinstall from the `git-workflow-autopush` memory if missing on a fresh
  clone). Never leave verified work uncommitted at a step boundary. Do not
  commit PDFs or `build/`.
- Figures: JFM canon (`post/closure_bound/figstyle_jfm.py` — STIX 8pt, boxed
  axes, inward ticks, no grid, italic panel labels, true 32pc width).
- Writing tone: plain, specific, measured; no self-congratulation, no AI-tell
  tics (see the global tone rule and `feedback-tone-no-self-congratulation`).
- Working style: pause after each step for a status + assessment; hold a
  reasoned coauthor view rather than deferring to every comment.
</content>
