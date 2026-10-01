# Changelog

## [0.2.2] - 2026-10-01

- With `symmetrize = true`, `factorized_yield` lets a mass yield whose complement is not
  tabulated stand for that complement as well, as it already did for `⟨TKE⟩(A)` and `σ_TKE(A)`.
  A yield measured on one wing only now gives both, as the docstring said; before, the unpaired
  mass was kept and its partner left without yield.
- `symmetrized_yield` likewise copies an unpaired cell `(A, TKE)` onto `(A₀ − A, TKE)` and lists
  the added masses. The sum grows by the yield of the unpaired cells; where every complement is
  measured it is unchanged.

## [0.2.1] - 2026-10-01

- `read_kinetic_energy_dispersion` accepts the optional uncertainty column,
  `A sigma_TKE sigma_TKE_uncertainty`, that a retrieval of `sigma_TKE_vs_A` writes. It refused
  such a file before.

## [0.2.0] - 2026-09-30

Breaking: a kinetic-energy width is no longer assumed.

- `factorized_yield` defaults to `dispersion = nothing`. A mass without a width is placed at its
  mean kinetic energy, shared between the two grid energies that bracket it so that the mean is
  exact; a mean outside the grid drops the mass. It was a Gaussian of 10 MeV.
- `read_kinetic_energy_dispersion` defaults to `default = nothing`, so a mass the table does not
  reach keeps no width. `kinetic_energy_dispersion` then returns `nothing`.
- `uniform_kinetic_energy_dispersion` requires its width.
- `DEFAULT_KINETIC_ENERGY_DISPERSION` is removed.

The width of `P(TKE|A)` belongs to the experiment that measured `⟨TKE⟩(A)`. Where that experiment
did not report one, the reconstruction now invents none.

## [0.1.9] - 2026-09-30

- `recommended_mean_total_kinetic_energy` and `RECOMMENDED_MEAN_TOTAL_KINETIC_ENERGY`: the energy
  standards of the mean pre-neutron TKE for ²⁵²Cf(sf) and thermal fission of ²³³U, ²³⁵U and
  ²³⁹Pu. They are Gönnenwein's 1991 recommendations as tabulated by Bertsch et al. (2015), with
  the ²⁵²Cf value from the absolute measurement of Henschel et al. (1981).

## [0.1.8] - 2026-09-30

- `factorized_yield(...; symmetrize = true)` imposes the exact pre-neutron identities:
  `Y(A) = Y(A₀ − A)`, `⟨TKE⟩(A) = ⟨TKE⟩(A₀ − A)` and `σ_TKE(A) = σ_TKE(A₀ − A)`. Where both
  complements are measured, each takes their mean; where one is, it serves both.
  `symmetrized_yield` does the same for a joint `Y(A, TKE)`. The default keeps the tables as
  measured.
- `factorized_yield` reads each mass once and sums in mass order. Unsymmetrized results agree with
  0.1.7 to rounding.

## [0.1.7] - 2026-09-30

- `read_segmented_curve`'s docstring gave `mean_of_ratios` as the averaging a per-charge partition
  inverts. It now names the charge-resolved relation, the only exact inverse, and the manifest's
  `[domain]` record. Documentation only.

## [0.1.6] - 2026-09-30

- The Gilbert–Cameron shell corrections ship with the package as
  `GILBERT_CAMERON_SHELL_CORRECTION_FILE`: Table III of the 1965 paper, pp. 1453–1455,
  transcribed three times independently. The transcriptions agree in every cell, and with
  Geant4's electronic copy in all 228 values. `S(Z)`
  past `Z = 98`, which the paper does not tabulate, is `NaN`. The reader treats a non-finite cell
  as untabulated, as it already treated trailing zero padding.
- Every test runs on a fresh clone. `FISSION_FRAGMENTS_DOMAIN_TEST_DATA` and the skipped-testset
  report are gone with the last test that needed them.

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
