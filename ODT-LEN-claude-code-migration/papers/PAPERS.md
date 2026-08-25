# Annotated bibliography — ODT & LEN

Substantial technical notes on the six project papers, written to be read directly by a
Claude Code session. Ordered by conceptual dependency, not date.

---

## 1. Kerstein, Ashurst, Wunsch & Nilsen (2001) — ODT vector formulation
*J. Fluid Mech.* **447**, 85–109. `full_text.pdf`

**Foundational model paper.** Generalises ODT from a single velocity component to a
**three-component vector velocity field** `v_i(y,t)` on a 1-D domain, which is what lets
the model represent pressure-scrambling (intercomponent energy transfer).

- **Two mechanisms:** molecular diffusion `(∂_t − ν∂_yy)v_i = 0`, `(∂_t − κ∂_yy)θ = 0`
  between events; a stochastic sequence of instantaneous **eddy events** for advection.
- **Eddy event:** `v_i(y) → v_i(f(y)) + c_i K(y)`, `θ(y) → θ(f(y))`.
  - **Triplet map** `f(y)`: shrink `[y0,y0+l]` to 1/3, three copies, middle reversed.
    Measure-preserving, continuous, scale-local. Chosen as the simplest map meeting these.
  - **Kernel** `K(y) = y − f(y)`: zero outside the eddy, integrates to zero (⇒ momentum
    conserved). Identity used: `∫K² dy = (4/27) l³`.
- **Pressure scrambling / return to isotropy:** amplitudes `c_i` chosen so total KE is
  conserved while KE is redistributed among components via a symmetric transfer matrix
  `T` (eq. 9) and a free parameter `α ∈ [0,1]`. `α=2/3` = equipartition of available
  energy; **`α=1` = maximum intercomponent transfer** and best matches DNS.
  "Available kinetic energy" `Q_i = (27/8) ρ₀ l v_{i,K}²`.
- **Eddy selection:** time scale from `(l/τ)² ∼ v²_{2,K} + α Σ_j T_{2j} v²_{j,K} − Z ν²/l²`
  (viscous penalty `Z`); rate `λ(y0,l;t) = C/(l²τ)`; sampled by a fixed-distribution
  **select-and-decide** (rejection) scheme so `λ` is never rebuilt.
- **Large-eddy anomaly:** rare huge eddies dominate transport (∝ size²); suppressed with
  a parameter-free **median-slope** method.
- **Validation:** temporally developing planar **mixing layer** and **wake** vs. DNS
  (Rogers & Moser; Moser et al.). Mean velocity, Reynolds shear stress, component
  variances, TKE budgets, and passive-scalar PDFs. Reynolds-stress/flux quantities are
  computed as **eddy-monitored fluxes**, not from `v_i` directly (Appendix).
- **Parameters used:** `Z=0.02`, `α=1` (case M) or `2/3` (case E); `C=3.78` (mixing
  layer), `5.55` (wake) — *note the 2005 erratum: correct C≈3.68.*

**Why it matters here:** this is the model your code implements. Every parameter (`C`,
`Z`, `α`) and every conserved-flux subtlety traces back to this paper.

---

## 2. Ashurst & Kerstein (2005) — variable-density ODT, mixing layers
*Phys. Fluids* **17**, 025107. `Ashurst_Kerstein___2005___...pdf`

Generalises ODT to **variable density** and introduces the **spatially developing
(S-flow)** formulation. ODT is described here as an outgrowth of the **linear-eddy model**
(Kerstein 1991) — relevant to the "LEN" question.

- **Density transport:** binary diffusion at fixed molecular number density gives
  `∂_t ρ = κ ∂_yy ρ` and `∂_t(ρ v_i) = μ ∂_yy v_i + κ ∂_y(v_i ∂_y ρ)`.
- **Second kernel `J = |K|`:** with variable density, `∫K dy = 0` no longer conserves
  momentum, so the eddy map becomes `v_i → v_i(f) + b_i J + c_i K` with `b_i = −H c_i`
  (H from density-weighted kernel integrals, eq. 22–23). `b_i = 0` recovers the
  constant-density case.
- **Available energy / selection** generalised (eqs. 24–29); same `α`, `Z`, `C` roles.
- **Compressibility (kinematic):** a **sonic-eddy** cutoff — disallow eddies with
  `Ma_eddy = l/(τ c_eddy) > M0` — to compare against compressible DNS at `Mc=0.7`.
  `M0=0.04`, `C=5.6` (temporal), `6.9` (spatial), `Z=0.02`, `α=1`.
- **Physics result:** entrainment of the **denser** stream is inhibited; layer growth
  should be **non-monotonic** in density ratio `s`, with a predicted trend reversal
  slightly beyond the experimentally explored range (`s≈8`). Compared to Brown & Roshko,
  Konrad, Pickett & Ghandhi experiments and Pantano & Sarkar DNS.
- **Appendix B:** the full **S-flow** control-volume derivation (mass/momentum/energy
  flux balances; the `V` lateral advection velocity; coordinate transform `y→ŷ`).

**Why it matters here:** the path to **variable-density jets and jet flames**, and the
formal spatially developing formulation used for round jets.

---

## 3. Lignell, Kerstein, Sun & Monson (2013) — adaptive mesh
*Theor. Comput. Fluid Dyn.* **27**, 273–295. `Lignell_et_al___2013___...pdf`

The **numerics** of the modern code: an adaptive, non-uniform, time-varying mesh with a
**Lagrangian finite-volume** advancement.

- **Triplet map on an adaptive mesh:** implemented exactly by displacing/duplicating cell
  faces (triples cell count in the eddy), avoiding the dispersion error of a fixed-grid
  `3ᵐ`-cell permutation (a size-6 permutation has only half the continuum mean-square
  displacement).
- **Mesh adaption:** distribute cells by **arc length** of the adapted profile(s)
  (parameter `Ndens` = points per unit arc length); constraints — min/max cell size and
  the **"2.5 rule"** (neighbour size ratio ≤ 2.5). Adaption happens after each eddy,
  after diffusion, and in long-untouched regions; it induces small numerical dissipation
  (KE loss ≈ 0.9% of physical at `Ndens=30`, 0.2% at 100).
- **Lagrangian advancement:** cells expand/contract with density (no mass crosses faces);
  continuity `ρΔx = const`; handles dilatation from heat release / spatial marching
  without an advective CFL limit. Full temporal and spatial conservation equations
  (continuity, species with Fick's law, `u/v/w` momentum, energy, a `dP/dt` pressure
  equation) are given in §5 — **useful as the reference for the code's discretisation.**
- **Demos:** turbulent channel (`Reτ=590`, `C=10`, `Z=600`; ~260 adaptive cells vs 38 M
  DNS cells → speedups of 8–68×), reacting ethylene temporal jet (`C=10, Z=200, β=0.9`),
  buoyant heated-wall spatial boundary layer (`C=3, Z=900, β=1`).
- **Appendix:** the **eddy sampling** (thinning/rejection), the `f(l)` size PDF, and the
  energy-based `1/τ = C√(2/(ρ₀l³))·(E_kin − E_pot − Z E_vp)` — the operational eddy-rate
  used in code.

**Why it matters here:** to modify or debug the C++ code you need this paper's mesh and
advancement model; the parameters and cost numbers set expectations for cluster runs.

---

## 4. Stephens & Lignell (2021) — the BYUignite/ODT C++ code
*SoftwareX* **13**, 100641. `full_text1.pdf`

The **software manual** for the reference implementation.

- **Repo:** `github.com/BYUignite/ODT`. MIT license. C++ (object-oriented) + Python.
  **Build:** CMake ≥ 3.12, **Cantera**, **yaml-cpp** (auto-built by the package), Git,
  optional Doxygen.
- **Architecture:** `main` → `domain` (owns objects, init via case-specific
  `domaincase`); three workers — **`solver`** (marches, invokes diffusion + eddies),
  **`micromixer`** (diffusive advancement: explicit Euler / semi-implicit / Strang,
  CVODE for stiff chemistry; talks to the `mesher`), **`eddy`** (samples/implements eddy
  events through rejection tests). Reacting cases use Cantera for properties/kinetics.
- **Directories:** `source/`, `build/`, `run/`, `input/<caseType>/input.yaml`,
  `input/gas_mechanisms/`, `data/` (raw + post), `docs/`.
- **Input:** human-readable `input.yaml`; parameters grouped (`params` class). Chemical
  mechanism in the input must match the one set at CMake configure time.
- **Running:** `run/` holds `odt.x` + scripts — `runOneRlz.sh` (1 realization),
  `runManyRlz.sh` (serial, set `nRlz`), and **SLURM** scripts `slrmJob.sh`,
  `slrmJob_array.sh` for parallel realizations (embarrassingly parallel; watch case
  names so realizations don't overwrite). Set `inputDir` and `caseName` in the run script.
- **Output & post:** `data/<caseName>/{input,runtime,data,post}`; Python tools in `post/`
  (per case type; `driver.py`), `post/tools/data_py.py` (text→binary numpy),
  `data_tools.py` (realization access, PDFs, domain bounds, etc.).
- **Example:** DLR-A canonical jet flame (`C=20, βLES=17, Z=400`, 1024 realizations).
- **Contributors named:** Sandia, **BTU Cottbus-Senftenberg** (Schmidt, Medina, Klein),
  Chalmers, CSE Inc. — i.e. the group this project sits in.

**Why it matters here:** this is your build/run/post workflow. When Claude Code edits
"the code," this is (a fork of) what it edits, and the SLURM scripts here are the
template for CARO/CARA jobs.

---

## 5. Medina Méndez, Sharma, Schmidt & Klein (2023) — far-field sound from ODT
*Proc. Appl. Math. Mech.* **23**:e202300186. `Proc Appl Math and Mech ... .pdf`

**The project's own north-star paper** (owner is a co-author). A framework to get
far-field **SPL** of a low-Mach jet from an ODT velocity field.

- **Analytical half:** start from compressible continuity/momentum/energy, form a
  **generalized nonlinear pressure wave equation** where pressure is an independent
  variable and the only nonlinearity is convective transport (Lighthill's spirit).
  Apply **low-Mach asymptotics** (`p* = p*₀ + Ma p*₁ + Ma² p*₂ + …`) with Strouhal
  `[S]=[Ma]⁻¹` (acoustic time scale). Orders: `O(Ma⁰)` ⇒ `p₀` constant; `O(Ma¹)` ⇒
  source-free wave equation for `p₁`; `O(Ma²)` ⇒ **wave equation with sources** for `p₂`
  (fluctuating Reynolds stress, turbulent heat flux, rate-of-change of KE / monopole,
  heat conduction). The `[S] ≠ [Ma]⁻¹` branch (turbulence time scale) is **deferred to
  future work**.
- **Spectra:** Fourier transform → Helmholtz eq.; Green's function + far-field
  approximation → spectral `p̂₁`, `p̂₂` **entirely determined by the velocity field**
  `V₀` (for isothermal jets, `T₀`, `ρ₀` constant). SPL from the pressure PSD
  `S_p = ⟨p̂† p̂⟩`, `SPL = 20 log₁₀(√S_p · p_b / p_ref)`.
- **Numerical half:** ODT supplies `V₀` for a **fully developed, statistically steady,
  locally homogeneous** round jet; the `∇·` operators reduce 3-D→1-D. Pressure isn't a
  variable in ODT — its scrambling is modelled by the kernel — and the paper argues the
  `p₁` effect is already captured by that kernel, so a **solenoidal (`∇·V₀=0`) vector ODT
  run suffices.**
- **Case (Table 1):** `Ma_j = 0.4`, `Re ≈ 553,000`, `Pr=0.71`, sideline 2.79 m, θ=165°;
  compared to Viswanathan (2007) experimental SPL and Tam LS/SS similarity spectra.
- **Result:** SPL agreement is good; `p₁` dominates most of the band, `p₂` dominates at
  high frequency. Framework extendable to variable-density flows and jet flames.

**Why it matters here:** this defines the deliverable pipeline (ODT velocity → spectral
sources → SPL) and names the exact validation case and its open ends.

---

## 6. Choi & Lumley (2001) — return to isotropy of homogeneous turbulence
*J. Fluid Mech.* **436**, 59–84. `thereturntoisotropyofhomogeneousturbulence 1.pdf`

The **physics that ODT's pressure scrambling emulates.** Hot-wire experiments on
grid turbulence after plane distortion, axisymmetric expansion, and axisymmetric
contraction, analysed with the **invariant technique**.

- **Return-to-isotropy tensor** `φ_ij` in the Reynolds-stress equation; general form
  `φ_ij = β b_ij + γ(b_ik b_kj + (2/3) II δ_ij)` (Cayley–Hamilton). **Rotta's linear
  model** = `γ=0`, `φ_ij = β b_ij`.
- **Key finding:** the return to isotropy is **nonlinear** — turbulence wants to become
  **axisymmetric** more than isotropic; trajectories in the `ξ–η` (invariant) plane are
  not straight lines to the origin. Rate is slowest for **cigar-shaped** (III>0)
  turbulence, fastest for **pancake-shaped** (III<0).
- **Turbulence triangle**, realizability conditions (Schumann), Reynolds-number
  dependence of the return rate `ρ*` — a full nonlinear model (eq. 58) is proposed.

**Why it matters here:** ODT applies "return to isotropy" **per eddy** through the kernel
`α`/transfer matrix — a much cruder device than Rotta/Choi–Lumley closures, which act on
averaged quantities. This paper is the yardstick for judging how faithfully the ODT
kernel reproduces real intercomponent transfer, and it connects the modelling to the
classical second-moment-closure literature (relevant when explaining or defending the
ODT pressure-scrambling choices).

---

## Add to the set
- **Sharma, Klein & Schmidt (2022)**, *Phys. Fluids* **34(8)** — far-downstream
  asymptotic velocity fluctuations in a round jet, an ODT study. The velocity-statistics
  basis for the acoustics work; referenced as [13] in the PAMM paper. Not currently in
  the Project — add it.
