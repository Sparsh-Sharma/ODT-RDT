# Project conventions

## Git commits
- Do NOT add `Co-Authored-By:` trailers to commit messages.

## Current work context — PRIORITY: paper-1 revision
The primary task on this branch (master) is fixing the rejected JFM paper
(JFM-2026-1781). Read `docs/REVISION_MAP_jfm1781.md` FIRST — it holds the
objection ledger, asset inventory, the combine-vs-split decision state, the
resubmission outline, and the ordered task list. Keep it updated as tasks
close or decisions land.

The Kerstein/Test-3 loop and the Gate A/B acoustics machinery live on
branch `LEN_Extension` (hand-off doc there: `docs/PLAN_paper2_test3.md`).
That work feeds the revision (see the map, section 3) but paper 1 outranks
it unless Sparsh says otherwise.

## Before any ODT development / simulation work — REQUIRED reading
Chats do not share history. The accumulated operational knowledge lives in
the project auto-memory (loaded as an index in every session started in this
directory). Before running, building, or campaigning ODT, READ these memory
files in full — do not rediscover their contents:
- `odt-rdt-project-state` — what the code variants are, branch/version map,
  which parameters exist (tStrainOn, allocMode, mapMidFrac, nSubKernelLevels),
  which campaigns ran and their verdicts, known bugs and their fixes.
- `hpc-cara-caro-usage` — where and how realizations run: caro via WSL ssh,
  build recipe (`build_caro/`), SLURM account/partitions, the node-packing
  script `run/caro/slrm_test3_node.sh`, and the quoting/stale-binary traps.
Campaign mechanics are also in `run/caro/README.md`. The LEN_Extension
worktree is at `..\ODT-RDT-len`; sessions started THERE get a different
(empty) auto-memory — start ODT-development chats from THIS directory.
