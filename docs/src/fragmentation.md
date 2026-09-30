```@meta
CurrentModule = FissionFragmentsDomain
```

# Charge distribution and domain

The isobaric charge distribution, the `Zₚ` model that supplies one where no evaluation exists, and
the set of fragmentations a model sweeps.

## Where a charge distribution comes from

`p(Z, A)` is a Gaussian in `Z` at each mass, so a model needs a centre `ΔZ(A)` and a width `σ_Z(A)`
over its whole mass range. Four things can supply them, in descending order of how much they are
worth, and [`build_charge_distribution`](@ref) resolves the last three in that order:

| | Source | Reached by |
|---|---|---|
| 0 | An evaluated table for this reaction | [`read_charge_distribution`](@ref) |
| 1 | The `Zₚ` model with this reaction's own fitted parameters | [`build_charge_distribution`](@ref), for the four reactions of [`WAHL_PER_REACTION`](@ref) |
| 2 | The `Zₚ` model with parameters estimated from the systematics | [`build_charge_distribution`](@ref), for any nucleus inside [`WAHL_VALIDITY`](@ref) |
| 3 | The conventional means, `ΔZ = −0.5` and `σ_Z = 0.6` | [`mean_charge_distribution`](@ref), or [`build_charge_distribution`](@ref) falling through |

Layer 0 is not part of the resolution: a table that is read is used, and the model is not
consulted. It is in the table because it is the thing the other three approximate, and because
a calculation that has one should use it.

The distinction between layers 1 and 2 is where the parameters of the same model come from. Table A
of the 1988 evaluation fits `ΔZ(140)`, `∂ΔZ/∂A'`, `σ_Z`, `F_Z` and `F_N` to the measured fractional
independent yields of one reaction; eq. (17) of LA-13928 [Wahl2002](@cite) estimates the same quantities from `Z_F`,
`A_F` and the excitation energy, for reactions that have no fit. Fig. 18 of that report puts the
cost of the estimate at a reduced `χ²` of 7.9 against 2.9. Measured against the four evaluated
evaluated tables, the fitted parameters halve the mean absolute difference in `ΔZ` for
²³⁹Pu and ²⁵²Cf and improve `σ_Z` for all four.

Which layer answered is recorded in the returned distribution's `source`, for a run to carry into
its metadata, because it is the difference between a result resting on this reaction's own data and one
resting on a trend across other reactions.

Energy is where the layers stop. All of this is fitted to spontaneous, thermal and resonance
fission; above roughly 8 MeV the 1980 and 2002 reports carry separate branches for the onset of
multi-chance fission, which this package does not implement and [`WAHL_VALIDITY`](@ref) excludes.

```@autodocs
Modules = [FissionFragmentsDomain]
Pages = ["charge.jl", "wahl.jl", "wahl_reactions.jl", "domain.jl"]
```
