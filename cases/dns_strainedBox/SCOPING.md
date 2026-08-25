# Strained-box DNS — scoping (caro)

Status: DRAFT for decision, 2026-08-25. Owner: Sparsh. Target cluster: **caro**
(via cara; SLURM `medium` partition, 256 CPUs/node, 42-day walltime cap).

## 1. Objectives (each maps to a concrete consumer)

| # | Objective | Consumer |
|---|---|---|
| O1 | Exact Π⁽ʳ⁾_ij = A_kl·M_ijkl from the full 3-D spectrum at several accumulated strains | Referee 1.1 (external benchmark); estimator "exact" reference |
| O2 | Azimuthally averaged Φ_mn → true a(κ₂,κ⊥), c(κ₂,κ⊥) **and** the discarded azimuthal harmonics D_mn | Direct A1-error measurement (manuscript §5.4's O(b²) claim) |
| O3 | Line spectra φ_mn(κ₂) + finite-ensemble LOS subsamples | `post/closure_bound/` estimator input; bound-width vs strain |
| O4 | Exact linear-RDT companion evolution of the same initial field, projected to 1-D | Referee 3.1 (correct projected-spectrum RDT reference) |
| O5 | HIT null case (S=0): estimator must report no strain signal | Alan's false-positive test; calibrates statistical noise floor |
| O6 | Polarization content γ_eff(e) = max|c|/a and departure of a from a(κ) | Adjudicates the stage-b hierarchy (polarization cap vs modulus structure) |

## 2. Flow configuration

Two-phase protocol (Lee & Reynolds 1985 lineage):

1. **Precursor**: forced HIT in a triply periodic box to statistical
   stationarity at target Re_λ; then forcing OFF and a short decay
   (~0.5 eddy turnover) to shed the forcing imprint.
2. **Strain phase**: impose homogeneous plane strain
   `A = S·diag(+1/2, −1/2, 0)` (stretch x₁ = chordwise, compress x₂ = the
   ODT line axis, neutral x₃), forcing off. Run to accumulated strain
   `e = S·t = 2.0` with spectral snapshots at
   `e ∈ {0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0}`.

Parameters:

- **Re_λ ≈ 100–140** at strain onset (enough inertial range for spectral
  claims; kinematics, not high-Re asymptotics, is the point).
- **Strain-rate sweep** S·k/ε at onset ∈ **{0.8, 4, 16}**: the paper's
  operating point, an intermediate, and a genuinely rapid case where linear
  RDT must be recovered (built-in verification against O4).
- **Ensemble**: ≥ 8 independent precursor realizations per S (different
  forcing seeds). Statistics: jackknife over realizations.
- **Null case**: same pipeline with S = 0 (O5).

## 3. Numerical method

Pseudo-spectral, Rogallo (1981) deforming frame. For irrotational plane
strain no remapping is needed within a run: wavevectors evolve as
`dκ_i/dt = −A_ii κ_i`, i.e. κ₁(t) = κ₁₀·e^{−e/2} (coarsening) and
κ₂(t) = κ₂₀·e^{+e/2} (refining); the momentum equation in the moving frame
gains only the linear term −A·u; 2/3 dealiasing; RK3/RK4 in time.

**Resolution management** (the classic strained-box issue): resolved κ₁
range shrinks by e^{e/2} (≈ 2.72 at e = 2). Strategy: **pad, don't remesh**
— anisotropic initial grid with N₁ = 3× base so the κ₁ coverage at e = 2
still exceeds the dissipation range; κ₂ coverage only improves. Remeshing
(with its interpolation losses) only becomes necessary for e > 2, which we
do not need.

**Small-scale resolution**: κ_max·η ≥ 1.5 at all times, checked per
snapshot (η shrinks under strain as energy is pumped in).

## 4. Grid & cost on caro

Base grid sized for Re_λ ≈ 120 precursor (κ_max·η ≈ 1.7 at 512³):

| Case | Grid (N₁×N₂×N₃) | Nodes (256c) | Precursor | Strain phase | Est. core-h |
|---|---|---|---|---|---|
| Pilot (Re_λ≈60) | 384×128×128 | 1 | ~2 h | minutes | ~1k |
| Production | 1536×512×512 | 4–8 | ~1–2 d | hours (S-dependent) | ~50–150k |
| Null (S=0) | 512³ | 2 | ~1 d | — | ~20k |

Notes: the strain phase is cheap (e = 2 at Sk/ε = 0.8 is ~2.5 eddy
turnovers; at Sk/ε = 16 it is ~0.1); cost is dominated by precursors, which
are shared across the S-sweep (same stationary states, different strain
launches — one precursor serves three S values + checkpointed restarts).
All well inside caro's `medium` limits. Storage: full spectral snapshots at
1536×512×512 ≈ 9.7 GB each (complex u_i, single precision post-hoc);
~7 snapshots × 3 S × 8 realizations ≈ 1.6 TB → keep on
`/gpfs/caro/scratch/ws/shar_sp-AssamTea` (DarjeelingTea holds ProLB);
reduce to (κ₂,κ⊥) objects (MBs) before transfer.

## 5. Code decision

No public pseudo-spectral code ships with plane-strain Rogallo deformation;
the modification is small (~200 lines in spectral space: moving wavenumbers,
−A·u term, time-dependent dealiasing mask). Candidates to modify:

| Candidate | Language | Pros | Cons |
|---|---|---|---|
| [hit3d](https://github.com/sthavishtha/turbulence-codes) (Stanford lineage) | F90 + FFTW + MPI | battle-tested HIT, GPL, cluster-native | slab decomposition (scales to ~N cores) |
| [HIT36](https://github.com/aroccon/HIT36) | modern Fortran | clean, actively maintained | smaller community |
| [spectralDNS](https://github.com/spectralDNS/spectralDNS) | Python/mpi4py/pyFFTW | fastest to modify, pencil decomp via mpi4py-fft | Python stack on caro to be verified |

**Recommendation**: prototype the strain terms in **spectralDNS** locally
(days, not weeks; validates the formulation against exact RDT at high S),
and in parallel check caro's module set (compiler, MPI, FFTW, Python) to
decide whether production runs use the Python code (if mpi4py-fft performs)
or a port of the validated terms into hit3d. Decision gate: pilot-case
throughput on 1 caro node.

## 6. Outputs specification (what leaves the cluster)

Per snapshot, computed in situ on the spectral field (small files only):

1. `phi_mn(k2)` — exact line spectra (κ₁,κ₃-integral; also = LOS-ensemble
   limit), all 6 components; plus subsampled finite-line ensembles
   (N_los ∈ {16, 64, 256}) for statistical error bars.
2. `a(k2,kperp)`, `c(k2,kperp)` — azimuthal averages (Batchelor–
   Chandrasekhar scalars), plus azimuthal harmonics m = 1, 2 of Φ_mn
   (the A1 error, O2).
3. `Pi_rapid_exact[3x3]` — direct quadrature of M_ijkl (O1).
4. `b_ij`, `k_t`, `ε`, `E(κ)`, κ_max·η.
5. RDT companion: same IC evolved under exact linear RDT (spectral ODEs,
   negligible cost), same outputs (O4).

Formats: HDF5 (or NetCDF), one file per (S, realization, e); naming
`sbox_S{S}_r{seed}_e{e}.h5`. Estimator ingestion: extend
`post/closure_bound/` with a reader mapping item 1 → `bound_pi` input and
items 2–3 → truth references.

## 7. Validation targets

- **Rapid limit** (Sk/ε = 16): b_ij(e) and component spectra vs exact RDT
  (O4) — must match to O(1/(Sk/ε)).
- **Lee & Reynolds (1985)**: b_ij evolution at matched S* (mind the
  factor-2 strain-parameter convention flagged by referee 1).
- **HIT limits**: precursor spectra vs standard forced-HIT references;
  estimator null test (O5).

## 8. caro reconnaissance (DONE 2026-08-25, via `ssh caro` from WSL)

- Login lands on `carologin6`; direct SSH works (keyed, `~/.ssh/id_rsa` in WSL).
- **Toolchain**: gcc 14.3.0, openmpi/4.1.7 (default, loaded), **fftw/3.3.10
  (loaded)**, hdf5/1.14.6, cmake/3.31.9, python/3.11.14 (loaded),
  miniforge3/25.3.0, intel-oneapi compilers + MKL available. Spack-user
  environment. Everything hit3d or spectralDNS needs is present.
- **Workspaces**: `AssamTea` and `DarjeelingTea` (GPFS) expire 2026-10-07
  (4 extensions left); `myfancyws` on /scratch expires 2026-10-09.
  → extend AssamTea before production; snapshots live there.
- **SLURM account**: `2002095` (used for all recent M03 jobs on `medium`);
  QOS available: short/long/interactive/vis/reservation.

## 9. Risks / open items

- [ ] Forcing-scheme imprint on anisotropy at strain onset → mitigated by
      decay gap; verify b_ij < 0.01 at onset.
- [ ] Single-precision snapshots sufficient for M_ijkl quadrature? (check
      on pilot against double.)
- [ ] Alan's view on the S-sweep values and on e_max = 2.
- [ ] AssamTea workspace extension before the production campaign.

## 10. Immediate next actions

1. ~~caro reconnaissance~~ DONE (§8).
2. Prototype spectralDNS + strain terms (locally or on carologin at small
   size), validate against exact RDT at Sk/ε = 64 on 128³.
3. Pilot on caro (1 node, 384×128×128, account 2002095, partition medium),
   decision gate on throughput.
4. Production sweep + estimator ingestion.
