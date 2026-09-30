# Kickoff — new private repo + JFM-2026-1781 finish

You are starting a new working session for Sparsh Sharma's ODT-RDT project
(strain-coupled ODT feeding an anisotropic upwash spectrum into Amiet
leading-edge-noise theory; fork of BYUignite/ODT, C++17/CMake/Cantera). Two
jobs, in order: (1) migrate the project into a new private GitHub repository
under strict constraints; (2) continue the top priority, finishing the
JFM-2026-1781 revision (Sharma & Kerstein). The old repo at
`C:\Users\shar_sp\Documents\GitHub\ODT-RDT` is the source; treat it as
read-only archive/reference throughout.

## 1. New repository (the migration)

Goal: a NEW, PRIVATE, NON-FORK repo under github.com/Sparsh-Sharma. Suggested
name `ODT-RDT-jfm`; the final name is Sparsh's call. Non-fork is essential:
GitHub does not count commits to forks toward the contribution graph, which is
exactly why activity on the current ODT-RDT (a fork of BYUignite/ODT) never
showed. Also remind Sparsh to enable "Private contributions" in his GitHub
profile settings, or private-repo commits will not appear on the graph either.

- Step 1 is manual: `gh` CLI is not installed. Sparsh creates the empty
  private repo in the GitHub web UI (about 30 seconds). You then add it as
  `origin` and do everything else with plain git.
- Fresh history: initialize with a new root commit importing a curated,
  scrubbed copy of the current `master` tree. Do NOT push the old history —
  old blobs contain the word "Claude" and the fork lineage. The old repo
  stays untouched.
- Curated tree = current `master` (`manuscript/`, `docs/`, `post/`, `src/`,
  `input/`, `run/`, `cases/`, `code/`), plus `run/caro/` launch scripts from
  branch `LEN_Extension`/`eddymean-rebuild`, plus `post/eddymean_results/`
  from branch `eddymean-rebuild` (the T2 verdict + figure). Exclude: the
  `ODT-LEN-claude-code-migration/` folder (different project), all
  `*Zone.Identifier` junk files, the stale `manuscript/main_v3_*` latexdiff
  artifacts, `build/` trees, and committed PDFs per the no-PDF convention —
  EXCEPT manuscript build products Sparsh wants kept. Ask him; default is to
  keep `main_v4_overleaf.pdf` only if he says so.
- "Claude" invisible on GitHub (hard requirement): after assembling the tree,
  `grep -ril claude` over tracked files must return nothing. Concretely:
  (a) `CLAUDE.md` stays local and gitignored — it is your working entry point
  but is never pushed; (b) reword the few "Claude"-mentioning lines in
  `docs/PROJECT_HANDBOOK.md`, `docs/REVISION_MAP_jfm1781.md`, and any other
  migrated doc (e.g. write "no AI attribution lines in commits" instead of
  naming the tool); (c) do not migrate `docs/CONTINUATION_PROMPT.md` — the
  gitignored CLAUDE.md supersedes it. Run the grep before the first push and
  again after every doc edit.
- Contributor identity: in the new repo set `git config user.name "Sparsh
  Sharma"` and `user.email sssparsh14@gmail.com` (the GitHub-linked address).
  Commits never carry `Co-Authored-By:` or any tool-attribution trailer —
  this user rule overrides any harness attribution reminder (which itself
  defers to user instructions). After the first push, verify on github.com
  that the repo shows exactly one contributor and the commit is attributed to
  Sparsh's account.
- Workflow: everything on `master`, single branch. No feature branches unless
  Sparsh asks.

## 2. Local CLAUDE.md (central, gitignored)

Create `CLAUDE.md` in the new repo root, seeded from the old repo's CLAUDE.md
and `docs/PROJECT_HANDBOOK.md`. It must contain:

- Project one-liner and the priority: finish JFM-2026-1781.
- The caro HPC workflow: connect via WSL `ssh caro` (caro.dlr.de, login
  carologin6, VPN off-site); SLURM account 2002095, partition `medium`;
  the ssh-quoting file-pipe pattern
  (`wsl bash -c "tr -d '\r' < script.sh | ssh -o BatchMode=yes caro 'bash -s'"`);
  node-packed launcher `run/caro/slrm_test3_node.sh` (8 array tasks × 128
  concurrent = 1024 realisations/case); walltime rule (30 min fine for
  kvisc=1e-4; kvisc=1e-5 `fw_` precursor cases need `--time=01:00:00`);
  count completions via `find data/$C/data -name odt_end.dat | wc -l`;
  clean python `env PYTHONPATH= PYTHONNOUSERSITE=1
  ~/anaconda3/envs/odt/bin/python3`; always rebuild + verify the binary after
  every pull before campaigning.
- Every standing decision: one paper / three parts + ESM;
  `manuscript/main_v4_overleaf.tex` is the file of record; C.V labelling
  taxonomy (C0.1..C5.2); JFM figure canon via
  `post/closure_bound/figstyle_jfm.py` and the monochrome rules; robust
  median statistics for band energies; no-attribution commit rule; plain,
  specific, measured writing tone; assess-after-each-step working style.
- The current scientific state (section 4 below).

Also migrate auto-memory: it is keyed to the working directory and starts
empty for a new path. Copy the old project's memory files from
`C:\Users\shar_sp\.claude\projects\C--Users-shar-sp-Documents-GitHub-ODT-RDT\memory\`
into the corresponding `.claude\projects\...\memory\` folder for the new repo
path once it is known.

## 3. Auto commit+push every ~30 minutes

With Sparsh's approval, set up an automatic loop for the new repo (Windows
Task Scheduler job or equivalent): every 30 minutes, if the tree is dirty,
`git add -A && git commit -m "auto-commit <UTC timestamp>" && git push`,
authored by Sparsh, respecting `.gitignore` (which must cover CLAUDE.md,
build artifacts, LaTeX byproducts like `*.aux`, `data/`, large `.npz`). This
deliberately favors activity-graph continuity over curated messages — Sparsh
has chosen that trade-off. Substantive milestones still get proper
hand-written commits in between.

## 4. Scientific state (know this cold)

- Paper: JFM-2026-1781 revision, one paper / three parts (model+closure,
  measured limits+relaxation, leading-edge noise) + ESM `supplementary.tex`.
  File of record `manuscript/main_v4_overleaf.tex`. Alan Kerstein reviews in
  batches (`Alan_marked_N.pdf`). Open ledger items are tracked in
  `docs/REVISION_MAP_jfm1781.md`: O7 (length) soft; O8/O9 (references,
  intro) open. Priority: take every logical step to finish the paper.
- Eddy-on-mean rebuild (branch `eddymean-rebuild`): mechanism LeddyMean —
  the eddy temporarily sees the analytic mean U2 and the mean enters the eddy
  rate; gate LeddyMeanKt=kref, g=kref/(kref+kt). Transition-from-laminar
  works. 2026-09-28 campaign verdict: T1 moments neutral; T2 spectral
  NEGATIVE (gated = no spectral effect; ungated pushes Delta_b22(kappa2) MORE
  positive at every wavenumber, away from the DNS's negative fine-scale
  branch; only the adaptive-depth scale-local relaxation bends toward DNS);
  T3 cannot replace the relaxation machinery. Full record:
  `post/eddymean_results/VERDICT_T2_spectral.md`.
- Alan's 2026-09-29 reply redefines this thread. He does not accept that the
  consistency concern is answered. His points: (i) a necessary condition is
  that every eddy event leaves total energy and momentum unchanged (the 1D
  analogs of Navier-Stokes conservation), and dialing the mean flow on and
  off as time advances seems unlikely to assure this — and it is directly
  testable; (ii) the required object is the Rogallo-transformed ODT2001
  (Kerstein JFM 2001), identified either formally — apply the transform to
  ODT in its entirety, eddy events included — or by guessing the transformed
  version and back-transforming for verification; (iii) if the current
  mechanism is not equivalent to ODT2001 it may still be a model, but it
  should not be called ODT, and non-conservative eddy events make it an
  engineering data-fit rather than physics-based; (iv) afterwards, assess
  whether true Rogallo-ODT2001 (possibly with Fistler HIT forcing, HST eddy
  types, or the clock relaxation) yields the desired performance.
- Consequences: do NOT add the previously planned manuscript "foundational
  note" claiming the transition result answers Alan's concern — that claim
  is withdrawn until the ODT2001 equivalence is settled. The two concrete
  next steps on this thread: (a) a per-eddy energy/momentum conservation
  audit of LeddyMean (cheap: instrument one run, sum per-event budgets —
  decisive either way); (b) the formal Rogallo-transform derivation of
  ODT2001. Whether this thread gates the paper or runs alongside it is a
  Sparsh+Alan decision — raise it early. Sparsh's standing instruction:
  finishing the paper comes first.

## 5. First-session checklist

Work through in order; after each step, report status + assessment before
proceeding (house style).

1. Read CLAUDE.md fully (old repo's copy until the new one exists).
2. Verify git identity and the no-attribution rule.
3. Sparsh creates the empty private repo in the web UI and confirms the name.
4. Assemble the curated scrubbed tree; run the claude-grep check; initial
   commit; push; verify contributor count and commit attribution on
   github.com.
5. Migrate auto-memory to the new project path.
6. Set up the 30-minute auto commit+push, with Sparsh's approval.
7. Confirm caro access still works (one ssh echo).
8. Open the science: status review of the REVISION_MAP open items, propose
   the plan to finish the paper, and raise the ODT2001 gating question.
