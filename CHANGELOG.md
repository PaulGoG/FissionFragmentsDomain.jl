# Changelog

## [0.1.5] - 2026-09-30

- `ManifestDomain` records whether the Gilbert–Cameron formula took its deformed branch,
  eq. (21) of the 1965 paper beside eq. (20), as `deformed_branch`. The field is required in a
  `[domain]` table with `level_density_model = "GC"`, and must be `false` or absent otherwise. The
  two settings give different `a` for about a third of the yield-weighted ²⁵²Cf fragments, and
  the record could not tell them apart.

## [0.1.4] - 2026-09-30

- The package version recorded in `ManifestDomain` is read once, when the package is compiled,
  instead of through `pkgversion` at every call. JET on Julia 1.11 reported a possible error
  inside that call, which failed the 0.1.3 test suite there. Behaviour is otherwise unchanged
  from 0.1.3.

## [0.1.3] - 2026-09-30

- `ChargeResolved`, the exact inverse of a partition that gives every fragmentation its own
  `a_L/a_H`: the relation is reduced over the charge distribution rather than through an
  effective ratio. `ChargeResolved(mean_total_excitation(masses, domain, ⟨TKE⟩))` also weights
  each fragmentation by its mean total excitation, which is the extraction's premise
  `ν ∝ E*` carried to mass-resolved multiplicities. `heavy_excitation_fraction`,
  `temperature_ratio` and `temperature_ratio_slope` take `(averaging, model, domain, A_H, ·)` for
  every averaging. The effective ratios miss this inverse by up to 10⁻² in `R_T` at
  `A_H ≈ 130`.
- The temperature-ratio run record gains an optional `[domain]` table, `ManifestDomain`. It holds
  the level density model, the ratio averaging and its excitation weighting, the charges per
  mass, the charge model, the mass table and the package version.
  `write_temperature_ratio_manifest` writes the record `read_temperature_ratio_manifest` reads.
- `charge_model_label`, `level_density_label`, `ratio_averaging_label` and `ratio_averaging` give
  the spellings for run records and configurations.
- Wahl (1988), CF252S: `ΔZ(A_F/2) = 0`, point X of Fig. 2. The steep branch reaches `ΔZ_max`
  only below `A_F/2`, and the model had carried `ΔZ ≈ 0.49` to symmetry, where a fragment is its
  own complement. No tabulated `A'` falls on the point, so Tables I–IV reproduce as before.

## [0.1.2] - 2026-09-30

- Every shipped test reference and fixture now rests on a public source.
  - The ²⁵²Cf fixtures take Wahl's 1988 model in place of a supplied charge table, so the domain
    and energetics tests run on a bare clone.
  - The archived-run yield projections, whose input table has no citable origin, are no longer
    shipped.
  - The comparison of the systematics with evaluated tables is now against the 1988 fits.
- `data/README.md` describes the charge distribution by its public route.

## [0.1.1] - 2026-09-30

- `charge_model` returns a `ChargeModel` on every path; where the systematics cannot form the
  precursor excitation energy it falls through to the means instead of returning `nothing`.

## [0.1.0] - 2026-09-30

First version: the fragmentation domain of `DeterministicSequentialEmission.jl`, extracted so that
the temperature-ratio extraction and the emission model share it.

- Nuclides, fissioning systems, the system token and its record; element symbols to `Z = 118`.
- Positional table reader; AME2020 mass excesses shipped, `Q`-values and separation energies.
- Isobaric charge distribution from evaluated tables, the Wahl (1988) model of four reactions, the
  Wahl (2002) systematics or the conventional means. The 1988 model reproduces that evaluation's
  Tables I–IV row by row from `ν̄_A` alone. Both `Zₚ` models use the lattice (erf) form of eq. (7)
  and serve fragments directly, with the proton even-odd factor and, on request, the neutron one; the fragmentation domain, with an optional vanishing polarization at the
  symmetric split, its charge windows and the split weight of a pair sum.
- Level density parameters, back-shifted Fermi gas and Gilbert–Cameron; their ratio across a pair
  averaged over the charge distribution in either order; the relation between the temperature
  ratio and the excitation-energy partition, with its inverse and slope.
- Energy balance and fragment kinematics.
- Fragment yields: joint `Y(A, TKE)`, its factorized reconstruction, the marginal `Y(A)`, and
  yield-weighted averages with uncertainty propagation.
- Curves of the heavy-fragment mass and the temperature-ratio run record, with the label of the
  systematic trend fixed as `systematic_trend`.
