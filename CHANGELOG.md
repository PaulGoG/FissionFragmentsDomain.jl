# Changelog

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
