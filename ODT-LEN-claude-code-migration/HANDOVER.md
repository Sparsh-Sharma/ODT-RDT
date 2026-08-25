# ODT & LEN — Project Handover

**Prepared:** 2026-08-23 · **Owner:** Dr. Sparsh Sharma (DLR Braunschweig)
**Purpose:** carry the full context of the *ODT & LEN* project from the claude.ai
Project into Claude Code, so a fresh session (human or Claude) can pick up without
re-reading everything.

---

## 1. What the project is

The project uses **One-Dimensional Turbulence (ODT)** — and the related
**Linear-Eddy-Model (LEM/"LEN")** lineage — as a **reduced-order, full-scale-resolving
stochastic turbulence model**, and couples it to an **acoustic analogy** to predict
**far-field sound radiation** from turbulent flows. The immediate target is
**low-Mach-number turbulent jets** (jet noise); the framework is meant to extend to
variable-density flows and jet flames.

The scientific bet: a DNS-quality turbulent velocity field is expensive in 3-D, but for
flows with a single dominant direction of statistical inhomogeneity (jets, wakes,
mixing layers, boundary layers) ODT reproduces the full range of turbulent length and
time scales on a **1-D domain** at a tiny fraction of DNS cost. That collapses the
acoustic source evaluation from a 3-D space→wavenumber + time→frequency transform to a
**1-D + time→frequency** transform (Medina Méndez et al. 2023). This is what makes
far-field spectra affordable.

This is a natural fit to the owner's research identity: stochastic turbulence modelling
(ODT) + LES/DNS/LBM + aeroacoustic theory, at DLR, in long-running collaboration with
the **BTU Cottbus-Senftenberg** ODT group (Heiko Schmidt, Marten Klein, Juan A. Medina
Méndez).

### Key prior results anchoring the project
- **Sharma, Klein & Schmidt (2022)**, *Physics of Fluids* **34(8)** — "Features of
  far-downstream asymptotic velocity fluctuations in a round jet: a one-dimensional
  turbulence study." (Owner's own ODT round-jet study; the velocity-statistics
  foundation for the acoustics work. **Not in the project PDF set — add it.**)
- **Medina Méndez, Sharma, Schmidt & Klein (2023)**, *PAMM* **23**:e202300186 —
  "Toward the use of a reduced-order and stochastic turbulence model for assessment of
  far-field sound radiation: low Mach number jet flows." (The acoustics + ODT
  framework; owner is co-author. **In the project set.**)

---

## 2. The model, precisely (so a new session can reason about it)

**Domain & advancement.** ODT evolves a three-component velocity vector `v_i(y,t)` (and
scalars: density/`ρ`, mixture fraction, enthalpy, species) on a 1-D line. Between eddy
events, only the truncated diffusion equation acts, e.g. `∂_t v_i = ν ∂_yy v_i`
(advective and pressure terms omitted). Two formulations:
- **Temporal** — the line evolves in time (homogeneous / temporally developing flows:
  channel, temporal jet, temporal mixing layer).
- **Spatial** — the line marches in a streamwise direction as a parabolic
  boundary-layer problem (spatially developing jets, wakes, plumes, mixing layers).
  Requires the advancement-direction velocity to stay positive everywhere.

**Eddy events = triplet map + kernel.** Turbulent advection is a stochastic sequence of
instantaneous **eddy events** on intervals `[y0, y0+l]`:
- **Triplet map** `f(y)`: compress the interval to 1/3, place three copies, invert the
  middle copy. It is measure-preserving (1-D analogue of solenoidal flow), continuous
  (no new discontinuities), and enforces scale locality; it conserves all
  domain-integrated moments (mass, momentum, energy).
- **Kernel operation** `v_i → v_i(f(y)) + c_i K(y)` (constant density) or
  `+ b_i J(y) + c_i K(y)` (variable density), with `K(y)=y−f(y)`, `J=|K|`. The
  coefficients redistribute kinetic energy among components subject to momentum and
  energy conservation, implementing **pressure scrambling / return-to-isotropy** at the
  level of the individual eddy (parameter `α`; `α=1` = maximum intercomponent transfer,
  the value that best matches DNS in the free-shear studies).

**Eddy selection.** Each candidate eddy `(y0,l)` gets a time scale `τ` from the
instantaneous **available kinetic energy** in the interval, minus a **viscous penalty**
(`Z ν²/l²`) that sets a threshold eddy Reynolds number. The rate distribution
`λ(y0,l;t) = C / (l² τ)` is sampled with a **thinning/rejection** ("select-and-decide")
method so the constantly-changing distribution never has to be rebuilt. Rare
unphysically large eddies are removed by a **large-eddy suppression** mechanism
(median method, scale-reduction method, or **elapsed-time** method `t > β τ`).

**Three tunable parameters:**
- `C` — eddy-rate coefficient → turbulence intensity / growth rate.
- `Z` — viscous-penalty coefficient → smallest allowed eddy (often irrelevant at high Re).
- `βLES` (or the suppression method's parameter) → largest allowed eddy.

**Numerics (modern code).** Adaptive, non-uniform mesh with a **Lagrangian
finite-volume** diffusive advancement (cells expand/contract; no mass crosses faces),
which also handles dilatation in spatial/variable-density flows without an advective
CFL limit. Mesh adaption merges/splits cells to keep an arc-length-based grid density,
subject to min/max cell size and a "2.5 rule" on neighbour size ratio. Stiff chemistry
via CVODE; Strang splitting available.

**Statistics.** Every observable is an **ensemble average over many independent
realizations** (each a different RNG seed). Reynolds stresses and budget fluxes are NOT
computed directly from `v_i` (the ODT velocities don't literally advect fluid); they are
monitored as eddy-induced fluxes so conservation laws hold exactly.

---

## 3. The literature set (the six project PDFs) and why each matters

Full annotations are in `papers/PAPERS.md`. In brief:

1. **Kerstein, Ashurst, Wunsch & Nilsen (2001), JFM 447, 85–109** — *ODT: vector
   formulation and application to free shear flows.* (`full_text.pdf`) The core
   three-component ODT formulation, pressure scrambling via the `K` kernel and the
   return-to-isotropy transfer matrix, eddy selection, the large-eddy anomaly. **This is
   the foundational model paper for the code you run.**
2. **Ashurst & Kerstein (2005), Phys. Fluids 17, 025107** — *Variable-density
   formulation and application to mixing layers.* Adds density transport, the `J` kernel
   for momentum conservation under variable density, the spatial (S-flow) formulation,
   and a sonic-eddy compressibility treatment. **The route to jet flames / variable
   density.**
3. **Lignell, Kerstein, Sun & Monson (2013), Theor. Comput. Fluid Dyn. 27, 273–295** —
   *Mesh adaption for efficient multiscale implementation of ODT.* The adaptive
   Lagrangian mesh and advancement scheme underpinning the modern code. **Read this to
   understand the code's data structures and numerics.**
4. **Stephens & Lignell (2021), SoftwareX 13, 100641** — *ODT: computationally efficient
   modeling and simulation.* The **BYUignite/ODT** C++ code paper: architecture
   (`solver`/`micromixer`/`eddy`/`domain`), `input.yaml`, run scripts, SLURM,
   post-processing. **This is your build/run manual.**
5. **Medina Méndez, Sharma, Schmidt & Klein (2023), PAMM 23:e202300186** — *Reduced-order
   stochastic model for far-field sound radiation: low-Mach jets.* The acoustics
   framework: low-Mach pressure wave equation, Mach-number asymptotics, SPL from the
   1-D velocity field. **This is the project's own north-star paper.**
6. **Choi & Lumley (2001), JFM 436, 59–84** — *The return to isotropy of homogeneous
   turbulence.* Experimental + invariant-technique treatment of return-to-isotropy
   (Rotta model, the turbulence triangle, realizability). **The physics ODT's pressure
   scrambling is emulating** — the conceptual check on the kernel's `α`/transfer model.

---

## 4. What is being migrated

Per the owner, the migration covers **more than the papers**:

- **ODT/LEN source code** — the C++ ODT solver (BYUignite/ODT and/or a BTU variant),
  plus any LEN/LEM code, plus Python post-processing and the acoustics/SPL pipeline.
  *Bring the git remote(s); do not copy build artifacts.*
- **Simulation data & results** — run outputs, cases, and post-processed data. These
  live on **CARO/CARA** and are large; they stay on the clusters. Claude Code works
  against **case definitions and post-processing scripts**, not the raw dumps.
- **Research notes / writeup** — derivations, drafts, and any paper-in-progress
  (e.g. the follow-on to the PAMM 2023 note). Bring these into `notes/`.

See `docs/recommended-structure.md` for exactly where each lands, and
`docs/claude-code-migration-guide.md` for how to connect the folder and the clusters.

---

## 5. Open threads / where the work stands (fill in / correct)

These are inferred from the literature and the stated scope — **confirm and edit**:

- The PAMM 2023 note is explicitly *"toward"* the framework and flags **further
  validation as future work**. Likely next steps: broader jet cases, the `[S] ≠ [Ma]⁻¹`
  Strouhal branch deferred in that paper, variable-density / heated jets, and directivity.
- Extending the acoustic-source evaluation beyond the isothermal constant-property
  (solenoidal) assumption toward variable-density jets and jet flames (needs the
  Ashurst–Kerstein 2005 machinery).
- Confirming/settling what **LEN** denotes in this project and whether it is a separate
  codebase to migrate.

---

## 6. First move in Claude Code

1. Set up the project folder and drop this kit into it (see
   `docs/recommended-structure.md`).
2. Put `CLAUDE.md` at the repo root.
3. Connect the folder in the Claude Code desktop client, then connect CARO and CARA
   (`docs/claude-code-migration-guide.md`).
4. Open a session and paste `prompts/first-session-prompt.md`.
