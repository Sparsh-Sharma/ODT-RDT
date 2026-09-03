# Digest: the two references Alan Kerstein pointed to (2026-08-26)

Both are ODT-native mechanisms for generating *anisotropy* inside the eddy
event, as alternatives to (or a principled home for) the explicit LRR/IP
rapid pressure–strain source the JFM manuscript bolts onto the standard
isotropic kernel redistribution.

## 1. Fistler, Kerstein, Wunsch & Oevermann, Phys. Rev. Fluids 5, 124303 (2020)
*Turbulence modulation in particle-laden stationary homogeneous shear
turbulence using ODT.* The single-phase HST case (Secs. II B–E, IV B–C) is
the relevant part; concept originally Scott Wunsch's.

**Setup.** ODT line along the mean-gradient direction x₂; initial uniform
shear; jump-periodic BC on the streamwise component, periodic on the others;
finite domain gives transition to statistical stationarity (domain size sets
the largest eddy and the stationary Re). No forcing needed — the mean-flow KE
reservoir is the source. Time normalized by shear collapses the transient.

**Deficiency of standard vector ODT (Sec. II E, verbatim in spirit):** the
component subject to production is distinguished, but "the other two
components ... are statistically equivalent" — "a partial representation of
HST anisotropy." It cannot capture "the structural difference between those
other components."

**The mechanism — three eddy types.** Each eddy event is idealized as a
*two-dimensional motion normal to one coordinate direction* (motivated by the
instability of uniform spanwise vorticity generating motion in the x₁–x₂
plane). Consequences:
- kernel energy exchange restricted to the *two in-plane components*
  (Eq. 12; the omitted component gets ΔE = 0);
- three eddy types indexed by the omitted component; sampling probabilities
  (Eq. 11; one parameter, 0.2, tuned to Rogers & Moin TKE partition) control
  the anisotropy;
- eddy time scale (Eq. 13) uses the available energy of the *two
  participating* components only;
- occurrence rates ≠ sampling probabilities in anisotropic flow, because the
  rate (Eq. 13) is state-dependent per type — a symmetry-breaking mechanism
  absent in isotropic flow, where equal probabilities give isotropy.
- "Three eddy types with equal sampling probabilities is a valid
  general-purpose alternative to ODT with one eddy type ... the unphysical
  statistical equivalence of two of the velocity components ... can then be
  broken by adjusting the probabilities without reformulating the model."

**Validation level:** Reynolds-stress partition vs Rogers–Moin (transient) and
Gualtieri et al. (stationary); Table I shows ODT b₃₃ overshoots (0.51 vs DNS
0.33). Empirical, as Alan says.

## 2. Kerstein, Fluids 7, 76 (2022), Sec. 6.2 "Allocation of Kinetic-Energy Changes"
Inside the AME (autonomous microscale evolution) framework. The section
generalizes the kernel energy redistribution into an *available-energy
allocation rule*:

- Component-i **available energy** Q_i = maximum KE extractable from
  component i by the kernels (minimize ΔE_i over c_i); Q = Σ Q_i.
- **Net available energy** H = Q + S_E, S_E any external energy
  source/sink. H ≤ 0 ⇒ eddy forbidden (criteria in Sec. 6.3).
- **Baseline ("no-memory") allocation:** final Q*_i = H/3 — "no basis for
  different outcomes for different components." Gives
  ΔE_i = S_E/3 + (1/3)(Q_j + Q_k − 2Q_i); for S_E = 0 this is the standard
  α = 2/3 redistribution. "A good baseline model and perhaps the only needed
  model."
- **Generalized allocation with memory parameter χ:**
  Q*_i = [1 + χ(3Q_i/Q − 1)] H/3,  χ ∈ [−½, 1];  χ = 0 recovers no-memory.
  Relation to the classic coefficient: χ = 1 − (3/2)α, α ∈ [0,1]; α is the
  fraction of Q_i transferred in equal shares to the other two components.
- **The point:** "regardless of the number and types of physical processes
  contributing to S_E, this input suffices ... to determine Q_i and H and
  thus ΔE_i" — a *unitary framework* in which every physical coupling enters
  the eddy event through S_E and one allocation rule.
- Momentum is *not* redistributed among components (vector constraint).

**Why "specific to isotropic flow":** the allocation is component-symmetric —
Q*_i depends on i only through Q_i. A mean-gradient tensor A_ij has no place
in it. Making the allocation direction-dependent (Q*_i informed by A_ij, or
by an eddy-type label as in Ref. 1) is the "more systematic approach to
shear turbulence and other anisotropic cases" Alan is gesturing at.

## What this means for the JFM reconstruction
- Our manuscript adds the rapid term as an explicit source *alongside* the
  standard α-redistribution and argues no double counting (Sec. 5.6). In the
  §6.2 language the rapid term is an S_E-like input and the *allocation* of
  H among components is where the strain tensor should act. That is a
  cleaner, ODT-native home for the closure than an external LRR/IP term.
- The DNS plane-strain result φ₁₁/φ₃₃ → −0.46 at e = 1 is precisely the
  u↔w structural difference that standard vector ODT cannot represent and
  that Ref. 1's three eddy types were built to capture. Our DNS quantifies
  the size of the effect Ref. 1 parametrizes empirically.
- Everything from the LOS estimator and the A1 control is *kinematic* — what
  a single line determines, and when the axisymmetric closure holds — and is
  independent of how ODT generates anisotropy. Those results constrain any
  ODT variant, including Alan's two routes.
