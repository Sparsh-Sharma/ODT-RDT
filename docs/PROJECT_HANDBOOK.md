# ODT-RDT — project handbook

Exhaustive reference for the ODT-RDT project (JFM-2026-1781 revision + the
eddy-on-mean rebuild). `CLAUDE.md` is the readable entry point and points here.
This file holds the full architecture, branch map, HPC recipes and gotchas,
results, and decisions. Verified against the repo on 2026-09-28; where a hash
or path is quoted, re-check with `git` before relying on it.

---

## 1. What the project is

Hybrid aeroacoustics: **Rapid Distortion Theory × One-Dimensional Turbulence →
Amiet leading-edge noise.** The frozen-isotropic upwash assumption in Amiet's
leading-edge-noise theory is replaced by a distorted, anisotropic upwash
spectrum produced by a strain-coupled ODT model (SC-ODT). ODT supplies the
line-resolved turbulence; RDT/strain distorts it; the distorted upwash feeds
Amiet's response to predict far-field sound.

- Code: fork of **BYUignite/ODT** (Stephens & Lignell 2021; C++17, CMake,
  Cantera). Remotes: `origin` = github.com/Sparsh-Sharma/ODT-RDT (push +
  fetch), `upstream` = github.com/BYUignite/ODT.
- Priority (Sparsh, blunt): **PAPER 1 — the JFM-2026-1781 revision — outranks
  everything.** Coauthor: Alan Kerstein (confirmed 2026-09-05). The paper is
  deliberately open-ended: new results land as paragraph upgrades, not
  restructures.
- LEN = Linear Eddy Model (context term from the sibling jet-noise project);
  "LEN_Extension" is a branch name, not the paper-2 label it sounds like.

---

## 2. Branch / worktree map (VERIFIED 2026-09-28)

`git worktree list` and `git branch -avv` are the source of truth.

| Ref | Tip | What it holds |
|---|---|---|
| `master` = `origin/master` (this worktree `ODT-RDT`) | a4d7dc6 | **THE PAPER.** `manuscript/` (main_v4_overleaf.tex + all figures + supplementary + rebuttals), `docs/REVISION_MAP_jfm1781.md`, `post/closure_bound/` analysis pipeline, strained-box DNS, `figstyle_jfm.py`. This working directory is on it — do paper work here. |
| `eddymean-rebuild` (origin + local) | 62dafbd | The eddy-on-mean **solver** + the T2 verdict (`post/eddymean_results/`), the campaign decks, `docs/HANDOFF_eddymean_campaign.md`, `docs/PLAN_rebuild_eddymean.md`. Branched off `LEN_Extension`. **All eddy-on-mean work belongs here — never on `master`.** |
| `..\ODT-RDT-pm` worktree | a4d7dc6 (detached) | origin/master materialised; additionally holds the six ~117 MB `fw_*_ensemble.npz` (gitignored) for re-running the eddy-on-mean spectral analysis. Keep for any deep-strain follow-up; otherwise redundant with `master`. |
| `..\ODT-RDT-len` worktree | 20ae1f9 `eddy-mean-forcing` | **A DIFFERENT active chat's parallel rebuild.** DO NOT TOUCH the worktree or the branch. |
| `LEN_Extension` (`origin/LEN_Extension` 76d2159) | — | Test-3 / clock-relaxation machinery; `run/caro/` launch scripts; `docs/PLAN_paper2_test3.md`. The caro clone tracks this branch. |
| `master-eddymean-backup` | 614a435 | Safety backup of the resolved 2026-09-28 master divergence (see below). Deletable once confirmed unneeded. |
| `Dev` (f8c1a01), `backup/master-pre-rewrite` (710de58) | — | Stale / pre-history-rewrite backup. |
| tag `jfm-2026-1781-submitted` | 4462f93 | Code state the rejected manuscript used (reproducibility anchor). |

**2026-09-28 master divergence — RESOLVED.** An earlier session committed
eddy-on-mean work directly onto `master` (an auto-commit locked it in at
614a435), so local `master` diverged from `origin/master` (ahead 85 / behind
158) and the `ODT-RDT` folder showed no `manuscript/`. Fixed by resetting
`master` to `origin/master` after backing the divergent line up to
`master-eddymean-backup`; the rebuild's unique content is also on
`eddymean-rebuild`, so nothing was lost. **Prevent recurrence: eddy-on-mean
work stays on `eddymean-rebuild`; `master` stays = the paper.**

**Solver-vs-analysis split (still true):** the rebuild SOLVER lives on
`eddymean-rebuild`; the paper's ANALYSIS pipeline + DNS live on `master`. The
spectral analysis runs the master-side `threeway.py` on ensembles produced by
the `eddymean-rebuild` binary — build the binary from `eddymean-rebuild` on
caro, run there, extract, then analyse in the `master` worktree (or
`..\ODT-RDT-pm`).

Manuscript file of record: **`manuscript/main_v4_overleaf.tex`** (on
origin/master). `main_v3.tex` retired 2026-09-16. Many `main_v3_*`/`main_v4_*`
diff and build artefacts also exist under `manuscript/` (latexdiff outputs vs
Alan's reviewed versions) — the live source is `main_v4_overleaf.tex`; ESM is
`supplementary.tex`.

---

## 3. HPC workflow — DLR caro cluster ("skills" not to re-derive)

All simulations run on **caro**, never locally. Runs before 2026-08-25 were
local Ubuntu/WSL; the decision to move everything to caro is standing.

### 3.1 Connection
- From Windows, drive it through WSL (not Windows-side ssh):
  `wsl bash -c "ssh -o BatchMode=yes caro '...'"`.
- `ssh caro` -> caro.dlr.de, login node **carologin6**, key `~/.ssh/id_rsa`.
  `ssh cara` -> cara.dlr.de, key `~/.ssh/id_ed25519`. (Historically caro was
  reached nested via cara/carologin3; direct WSL access has worked since
  2026-08-25 — use it.) Needs Cisco AnyConnect VPN when off-site. Job mail:
  sparsh.sharma@dlr.de.
- SLURM account **2002095**. Partition **`medium`** (default, 256 CPUs/node,
  257 GB, whole nodes, `--exclusive` fine). Idle short partitions `ppp`
  (1 node, cpu=32/mem=30G QOS cap — never `--exclusive`) and `grace` (ARM Grace,
  the x86 conda `odt.x` will NOT run there — avoid) are sometimes free when
  `medium` is queued; check `sinfo` and `squeue -j ID -o "%S %R"`.

### 3.2 The SSH-quoting trap (bit this repeatedly)
Shell variables and `for` loops inside inline `ssh caro '... $VAR ...'` arrive
EMPTY or mangled on the remote side; quoted heredocs break too. **Reliable
pattern:** write the script to a file, then pipe it over stdin:
```
wsl bash -c "tr -d '\r' < /path/script.sh | ssh -o BatchMode=yes caro 'bash -s'"
```
Pass only literal values in inline `sbatch` lines. `scp` from `/mnt/c` carries
CRLF — strip with `tr -d '\r'` (or `sed 's/\r$//'`) before dropping into the
Linux tree. Never pipe `git pull` through `tail`/`grep` (it hides "Aborting" on
modified/untracked files); echo `PULL_EXIT=$?` and verify `git log -1`.

### 3.3 Build on caro
Clone/worktree at `~/ODT-RDT` (and worktrees like `~/ODT-RDT-em` for
eddymean-rebuild). CMake against the pre-existing conda env
`~/anaconda3/envs/odt` (Cantera + yaml-cpp + conda gcc; RPATH baked in, no
LD_LIBRARY_PATH):
```
git worktree add ~/ODT-RDT-em eddymean-rebuild
cmake -S . -B build_em -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=$HOME/anaconda3/envs/odt \
  -DCMAKE_CXX_COMPILER=$HOME/anaconda3/envs/odt/bin/x86_64-conda-linux-gnu-c++
cmake --build build_em -j8
cp build_em/src/odt.x run/odt.x
```
(Or `bash build_odt.sh` in the worktree.) **ALWAYS rebuild + verify (e.g. posf
at first dump) after pulling code, before campaigning** — a caro script once
pulled but did not rebuild and a whole campaign ran a stale pre-feature binary.
`odt.x` exits(0) silently if the target `data/CASE/data/data_NNNNN` already
exists — do not pre-run smoke tests into the campaign output dir.

### 3.4 Node-packed launcher and case size
`run/caro/slrm_test3_node.sh` (on LEN_Extension / eddymean-rebuild): 8 array
tasks per case, each running 128 concurrent `odt.x` on one node via
`xargs -P 128`. A case = 8 nodes × 128 = **1024 realisations**. Submit:
```
sbatch --account=2002095 --partition=medium --export=CASE=<case> slrm_test3_node.sh
```
`medium` allocates WHOLE 256-core nodes regardless of `--cpus-per-task`, so
never submit 1-realisation-per-task arrays there (256× the node-hours). Seeds
are paired across cases (`seed = 22 + shift`, shift 0..1023) so contrasts are
same-seed. `run/caro/slrm_test3_array.sh` is the older naive array (avoid on
medium).

### 3.5 Walltime sizing (important, learned this cycle)
- The launcher default `--time=00:30:00` is correct for kvisc=1e-4 cases
  (e.g. `homogeneousStrain2`, ~3 min packed).
- The **kvisc=1e-5 `fw_` precursor cases** (tStrainOn=0.4, gate-A/allocation
  protocol) are BACK-LOADED: they crawl through the precursor transient to
  sim-time t~0.09 (~15-18 min, ~24M eddy trials) then accelerate through the
  strained phase. Solo/exclusive `fw_S40I` (tEnd 0.425) ~23 min; packed ~35 min.
  **Resubmit `fw_` cases with `--time=01:00:00`** — the 30-min default KILLS
  them mid-transient with zero completions (this happened and cost a full
  48-node launch).
- To size an unknown case, run one solo probe:
  `sbatch --nodes=1 --exclusive --time=01:30:00` a single
  `./odt.x <case>_probe 0` under `/usr/bin/time -v`.

### 3.6 Monitoring
- True completions: `find data/$C/data -name odt_end.dat | wc -l`. NOT
  `ls data/$C/data | wc -l` — that counts STARTED realisation dirs (= 1024 from
  launch, misleading).
- Queue: `squeue -u shar_sp`. After any SLURM controller outage, re-verify job
  ownership (`scontrol show job` UserId) — controller resets have re-issued job
  IDs to other users; dead jobs leave frozen runtime files with no error.
- `S40`/rapid cases: a known ~0.4% `domainPositionToIndex` abort at high strain
  (eddy position ≈ ±4 L from dilatation). Harmless — the extractor NaN-masks it
  and the analysis uses medians. `/tmp` is login-node-local (different
  carologinN per session) — do not leave logs there.

### 3.7 Extract → analyse
- Extract ensembles: `post/extract_ensemble.py` with env
  `ODT_WS=<worktree> ODT_CASE=<case>`, run with the clean python of §3.3. Output
  `<ws>/<case_lower>_ensemble.npz` (keys `lines` (n_dumps, n_rlz, 2048, 3),
  `Ldump`, `times`). (`post/closure_bound/odt_null/extract_ensemble.py` is a
  sibling copy.) Transmission bands:
  `post/acoustics/results_test3/dump_bands.py` -> `bands_<case>.npz`.
- Runs read input from `../data/<case>/input/input.yaml` (NOT `../input/`);
  dumps to `../data/<case>/data/data_NNNNN/dmp_*.dat`.
- scp the npz home to the `..\ODT-RDT-pm` analysis worktree, then delete the
  caro copy (next `git pull` there otherwise collides with the committed copy).

---

## 4. Analysis pipeline (on origin/master / `..\ODT-RDT-pm`)

`post/closure_bound/` holds the paper's analysis. Key scripts:
- `strained/threeway.py` — `observables()` bins Δb_nn(κ2) in κ2(e)/κc(0), then
  compares ODT / exact-projected-RDT / DNS. `dns_like()` reads the 128³ DNS
  from `strained/n128/chk_r{0.8,16}_s*_e{0,1}.npz` (keys `kappa2`, `phi_line`).
  Produces `fig_threeway`.
- `strained/fourway.py` (+ `fourway_rcs.py`, `refresh_threeway.py`) — consumes
  `fw_*_ensemble.npz`; adds the ODT + adaptive-depth scale-local relaxation
  (ODTK) row. `fourway.npz` / `threeway.npz` are archived reductions
  (`threeway_archived.npz` = pre-dilatation-fix ODT rows). Produces
  `fig_fourway`.
- `robust_optionA.py` — median-of-ratios + paired-seed contrasts for the band
  npz. **Band energies are heavy-tailed — always medians, never means**
  (single realisations sit 100-2000× above the median; means never converge
  even at 1024 rlz).
- `figstyle_jfm.py` — the JFM figure canon (STIX 8pt, boxed axes, inward ticks,
  no grid, italic panel labels, built at true 32pc = 5.31 in, never scaled).
- Wavenumber units differ: ODT line κ is physical (up to 1e3+); DNS κ2 is
  integer modes (1..64). Centroid-normalise both (`centroid()`) — a raw overlay
  pins DNS b22 at −1/3 (garbage).
- Rapidity matching: match ODT χ = S·k_t/ε to DNS 0.8 / 16 via the S05/S8 (a.k.a.
  S2/S40 precursor) configs; `homogeneousStrain2` alone is χ~1.2.

The eddy-on-mean T2 driver written this cycle is `fw_eddymean.py` (in
`ODT-RDT-pm/post/closure_bound/strained/`; a copy is committed under
`post/eddymean_results/` on eddymean-rebuild). It calls the threeway
`observables()` on the fw ensembles.

---

## 5. Paper state — JFM-2026-1781 (Sharma & Kerstein)

Read `docs/REVISION_MAP_jfm1781.md` on origin/master for the authoritative
ledger. Summary:

- Rejected 18 Aug 2026, 3 referees. Resubmission of substantially the same
  paper to JFM is barred; any future JFM submission with original material must
  cite JFM-2026-1781 and returns to the same associate editor (Sutanu Sarkar).
  So the paper must be visibly reconstructed and the cover letter must map each
  objection to a change. All three referees found the §5 closure the
  interesting part — it is the spine.
- **Architecture (decided 2026-09-08, per editor guidance): ONE paper, THREE
  parts + ESM.** Part I = the model and what a line can know (formulation +
  closure/bounds + measured A1 + three-way DNS validation). Part II = measured
  limits and the scale-local relaxation mechanism (transmission diagnostic,
  interventions ledger, the four/five-way Δb22 figure). Part III = leading-edge
  noise consequence (Amiet mapping from the fitted RDT-vK family, ΔSPL vs strain
  and Sk/ε with the validity envelope). ESM = CMK validation, level1a numerics,
  full intervention family.
- Objection ledger status (O1–O9): O1/O2 evidence done (three-way DNS + exact
  RDT; rigid-translation claim corrected — the eddy-free run is ODT kinematics,
  exact projected RDT is the linear reference); O3 reframed to exact-obstruction
  + LP bounds + measured A1 + exact-for-axisymmetric; O4 partial; O5 diagnosed
  (component ordering lives in the LRR b-term, not ODT); O6 = the leading-edge
  section (D1 was settled to Option S then re-folded into the one-paper Part
  III); O7 length (main ~41 pp + 5 pp ESM, target was ~30 — soft), O8/O9
  references/intro OPEN. Check the map for the live status.
- Manuscript file of record: `manuscript/main_v4_overleaf.tex` (the promoted
  compact version, ~47pp/8pp at promotion 2026-09-18; ~41 pp main after Part
  III). Terminology: category.variant C0.1..C5.2 taxonomy (headline names:
  standard ODT / SC-ODT / clock-relaxed ODT); SC-ODT reference pair = C0.
- Alan proofreads in batches (`Alan_marked_N.pdf`); implement + latexdiff vs his
  reviewed version. Manuscript diffs live under `manuscript/main_v4_alan*_diff.tex`.

Key quantitative results already in the paper (from the memory + notes):
- 128³ strained-box DNS (Rogallo deforming frame, decay precursor, e≤1, 4
  seeds, Sk/ε ∈ {0.8, 16}) is the production benchmark; 256³ shelved (too
  expensive for one 24h node). Single-line closure EXACT for axisymmetric
  distortion; its plane-strain A1 error is a measured strain-geometry property.
- b22(e) (the acoustically relevant upwash amplification) is linear to e=1 and
  geometry-robust.
- RDT-distance envelope: rms distance of strained ODT line spectra from exact
  Cauchy-RDT-distorted vK grows 13% at onset → ~22% (Sk/ε≈0.4) vs ~24%
  (Sk/ε≈8) at e=2 — departure grows with strain rapidity, localised in the
  rapid kinematics (κ-uniform operator + rigid dilatation).
- Four-way Δb22(κ2) at e=1: baseline SC-ODT is flat/positive (the
  wavenumber-uniform deficiency Referees 1/3 flagged); the adaptive-depth
  scale-local relaxation (ODTK) tilts toward the DNS (removes ~half the
  fine-band misallocation at Sk/ε=0.8); at Sk/ε=16 the DNS lies on exact RDT and
  the relaxation moves it only marginally within e≤1 (resupply saturation).
- Part III: distorted inflow lowers far-field SPL by ~5/9/14 dB at e=0.5/1/2
  across the resolved band; four-way experimental validation against Paterson,
  Bampanis, Narayanan cases (chain in `post/acoustics/expval/chain.py`).

---

## 6. Eddy-on-mean rebuild (this cycle's main result)

Answers Alan Kerstein's objection (2026-09): a proper turbulence model must
generate turbulence from a laminar (u'=0) start under strain alone — the mean
gradient must enter the eddy machinery (his Rogallo-frame consistency argument),
not only the continuous production term.

**Mechanism `LeddyMean`** (default off = bit-identical): each eddy attempt
temporarily adds the analytic mean U_2 = A_22 (y − xc) to `vvel`, runs the
normal eddy (rate `eddyTau`, triplet map, kernels see the full velocity), then
removes it — leaving the sawtooth (M−I)[U_2] and putting the mean into the eddy
RATE. One eddy of size l injects (4/27)(a l)² into the upwash. R-independent, so
it SUPPLEMENTS the exact production P_22 = 2 a R_22 (kept in FULL), not replaces
it: at k_t→0 only the eddy fires (transition); once turbulent, production
dominates. Code touch: `src/solver.{h,cc}` (wrapper `shiftLineMean`, `meanGate`),
`src/domain.cc` (laminar guards for LRR 0/0 and singular Lyapunov at R=0;
production kept), `src/param.*`,
`src/domaincases/domaincase_odt_homogeneousStrain.cc` (LlaminarIC). Params:
`LeddyMean` (bool, off), `LlaminarIC` (bool, off; u=v=w=0 transition test),
`LeddyMeanKt` (double, gate scale, default 1e-3).

**Gate `LeddyMeanKt` = kref:** injected mean scaled by g = kref/(kref+kt),
k_t = line TKE. g→1 laminar (transition fires); g→0 once turbulent (reverts to
normal ODT, no over-isotropising). Without it the mean over-fires eddies and
degrades b22 (0.095→0.049 vs DNS 0.104 at c=1.45); gated recovers to ~0.076.
OPEN issue: kref is a NEW constant, against the paper's "no new constants"
ethos. Parameter-free target (sequel work): fire a mean-driven eddy only when
the fluctuation alone could not afford one (compare eddyTau energy with/without
mean).

**Campaign 2026-09-28 (commit 62dafbd; binary 2f7838a):** 6 × 1024-rlz `fw`
ensembles — gated EM + ungated-probe EMU (LeddyMeanKt=1e9) + same-binary
baseline I — at Sk/ε~0.8 (`fw_S2`) and ~16 (`fw_S40`), precursor protocol,
kvisc 1e-5. Δb22(κ2) at e=1 vs 128³ DNS + exact RDT via the threeway pipeline
(`fw_eddymean.py`).

**VERDICT (`post/eddymean_results/VERDICT_T2_spectral.md` on eddymean-rebuild):**
- **T1 moments: NEUTRAL** — gated recovers baseline/DNS b22, no degradation.
- **T2 spectral: NEGATIVE (decisive).** DNS/RDT target bends negative at fine
  scales (upwash below the transverse mean — the branch the referee-1/3
  wavenumber-uniform objection needs). Baseline is flat +0.08 (the deficiency).
  Gated EM sits ON the baseline (no effect — gate off once turbulent, by
  design). Ungated EMU pushes b22 MORE positive at every wavenumber (+0.10 at
  Sk/ε=0.8, up to +0.29 at 16) — the WRONG way, away from the DNS. Only the
  archived adaptive-depth relaxation (ODTK) bends the fine scales toward the DNS.
- **T3 interventions: eddy-on-mean CANNOT replace the adaptive-depth/clock
  relaxation** (wrong sign), so it does not simplify the paper.

Per-band Δb22(κ2) at e=1 (band edges κ2(e)/κc(0) 0.30..12):

```
Sk/eps ~ 0.8:
  DNS   +nan +0.032 -0.033 -0.029 -0.032 -0.036 -0.040 -0.046   (bends NEGATIVE)
  I     +0.089 +0.103 +0.103 +0.096 +0.073 +0.075 +0.079 +0.076 (flat positive)
  EM    +0.105 +0.103 +0.066 +0.084 +0.073 +0.070 +0.067 +0.080 (= baseline)
  EMU   +0.147 +0.135 +0.120 +0.118 +0.100 +0.095 +0.086 +0.097 (wrong way)
  ODTK  +0.076 +0.094 +0.061 +0.054 +0.037 +0.040 +0.037 +0.024 (toward DNS)
Sk/eps ~ 16:
  DNS   +nan +0.143 +0.007 -0.014 -0.037 -0.046 -0.049 -0.053
  I     +0.124 +0.114 +0.123 +0.115 +0.113 +0.111 +0.109 +0.107
  EM    +0.127 +0.117 +0.121 +0.120 +0.120 +0.115 +0.114 +0.116
  EMU   +0.288 +0.254 +0.246 +0.251 +0.226 +0.203 +0.177 +0.158
  ODTK  +0.123 +0.119 +0.115 +0.110 +0.108 +0.096 +0.097 +0.099
```

**DECISION (recommended, pending Alan):** keep the paper's spectral results
as-is (they rest on scale-local relaxation, which works); fold eddy-on-mean in
only as a short FOUNDATIONAL note — the formulation now transitions from laminar
and embeds the mean in the eddy, answering Alan's consistency concern — WITHOUT
claiming a spectral gain. The parameter-free version is sequel work. Note to
Alan drafted at `notes/email_alan_eddymean.html` (untracked). Figure
`post/eddymean_results/fig_fw_eddymean.{pdf,png}`. Value of the rebuild:
transition capability at spectral neutrality, i.e. foundational correctness, not
a better result.

Two rebuild/campaign hand-off docs live on the `eddymean-rebuild` branch (read
via `git show eddymean-rebuild:docs/...` or in a worktree of that branch):
`docs/HANDOFF_eddymean_campaign.md` (the three-test brief T1/T2/T3, code state,
pitfalls) and `docs/PLAN_rebuild_eddymean.md`.

---

## 7. Conventions and gotchas

- **Commits: NO `Co-Authored-By:` / no "Generated with Claude" line — ever.**
  Author is Sparsh Sharma alone (hard rule; history was rewritten to enforce
  it, 2026-09-03). This overrides any harness default attribution guidance.
- A local `.git/hooks/post-commit` auto-pushes the current branch after every
  commit (also catches GitHub Desktop commits; non-fatal offline;
  `push.autoSetupRemote=true`). Reinstall from the `git-workflow-autopush`
  memory on a fresh clone. Never leave verified work uncommitted at a step
  boundary. Do not commit PDFs or `build/`. `backup/master-pre-rewrite` holds
  the pre-rewrite history.
- Figures: JFM canon via `post/closure_bound/figstyle_jfm.py`; build at true
  \textwidth, never scale down; consistent colour palette across the paper.
- Writing tone (hard rule): plain, specific, measured. No self-congratulation,
  no rule-of-three triads, no "not just X but Y", no empty intensifiers, no hype
  openers/closers. Let the science stand.
- Working style: assess after each step (status + assessment before the next
  action). Hold a reasoned coauthor view; recommend rather than ping-pong to
  every comment.
- Environment: heredocs with backslashes mangle in this Windows/PowerShell +
  WSL setup — write edit scripts via the Write tool, not inline heredocs. TinyTeX
  / latexmk compiles the manuscript locally.
- Untracked-by-design: `notes/*` (email drafts to Alan, e.g.
  `email_alan_eddymean.html`) and `JFM_2026_ODT_RDT.zip` (manuscript lineage)
  are intentionally not committed.

---

## 8. Immediate next options (as of 2026-09-28)

1. Send Alan the eddy-on-mean results note (`notes/email_alan_eddymean.html`);
   get his read on folding the transition result in as a foundational note.
2. On Alan's verdict, draft the short foundational-note paragraph into
   `manuscript/main_v4_overleaf.tex` (on origin/master / `..\ODT-RDT-pm`) — no
   spectral claim.
3. Continue the open JFM-1781 ledger items (O7 length, O8/O9 references/intro)
   per `docs/REVISION_MAP_jfm1781.md`.
4. Sequel work: the parameter-free eddy-on-mean (fire a mean-driven eddy only
   when the fluctuation alone could not afford one) — removes the `kref` constant.
5. Housekeeping: `master-eddymean-backup` can be deleted once Sparsh confirms the
   old divergent auto-commit is not needed (its content is on `eddymean-rebuild`).
