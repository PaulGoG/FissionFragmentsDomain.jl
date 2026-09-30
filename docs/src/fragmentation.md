```@meta
CurrentModule = FissionFragmentsDomain
```

# Charge distribution and domain

The isobaric charge distribution of primary fragments, the `Zₚ` model that supplies it, and the
set of fragmentations a model sweeps.

## Where a charge distribution comes from

A [`ChargeModel`](@ref) gives `p(Z | A)` for fragments of every mass. Four sources can supply one,
in descending order of how much they are worth, and [`charge_model`](@ref) resolves the last
three in that order:

| | Source | Reached by |
|---|---|---|
| 0 | An evaluated table for this reaction | [`read_charge_distribution`](@ref) |
| 1 | The `Zₚ` model with this reaction's own least-squares parameters, Wahl (1988) [Wahl1988](@cite) | [`Wahl1988`](@ref), for the four reactions of [`WAHL_1988`](@ref) |
| 2 | The `Zₚ` model with parameters from the systematics across reactions, LA-13928 [Wahl2002](@cite) | [`WahlSystematics`](@ref), for any nucleus inside [`WAHL_VALIDITY`](@ref) |
| 3 | The conventional means, `ΔZ = −0.5` and `σ_Z = 0.6` | [`mean_charge_distribution`](@ref) |

Layer 0 is not part of the resolution: a table that is read is used, and the model is not
consulted.

The distinction between layers 1 and 2 is where the parameters of the same model come from. Table
A of the 1988 evaluation fits them to the measured fractional independent yields of one reaction.
Eq. (17) of LA-13928 estimates them from `Z_F`, `A_F` and the excitation energy, for reactions
that have no fit of their own. Fig. 18 of that report puts the cost of the estimate at a reduced
`χ²` of 7.9 against 2.9.

## The 1988 model, reproduced

[`Wahl1988`](@ref) is the 1988 evaluation's own calculation, not an approximation to it. Given
only the tabulated `ν̄_A`, it returns every `Zₚ(A)` of the evaluation's Tables I–IV to within
0.0007, and every calculated `Z̄(A)` and `RMS(A)` to within 0.0018, for all four reactions.
[`fractional_independent_yields`](@ref) is the product-level form those tables tabulate.

Two readings of the paper are needed for that, and both are settled by the tables themselves,
not by assumption:

- **The near-symmetry geometry of Fig. 2.** The steep branch of `ΔZ` crosses `ΔZ = 0` on the
  `Zₚ = 50` line and is displaced by `ΔA'_Z` towards higher `A'` at the height of the peak line.
- **The regions.** `σ̄₅₀` applies on the steep branch and `F = 1` below the junction, both in
  `A'`. Footnote a of Table A lists the same regions as product masses, and misprints one: the
  U235T range ends at 108, not 109.

The parameters are those the tables were calculated with, each within the rounding of the
printed value; [`WAHL_1988_TABLE_A`](@ref) keeps the values as printed.

## Fragments, not products

Wahl's yields are yields of products, after prompt-neutron emission. The model is a function of
the precursor mass `A'`, so it applies to primary fragments at `A' = A` directly. One thing does
not carry over: the neutron even-odd factor. It was fitted to the parity of the product's `N`,
and a fragment's parity differs from its product's by the number of neutrons it evaporates.
[`fragment_charge_yields`](@ref) therefore applies the proton factor alone, unless the model is
built with `neutron_pairing = true`.

## Energy

Energy is where the layers stop. Everything here is fitted to spontaneous, thermal and resonance
fission. Above about 8 MeV, the 2002 report carries separate branches for the onset of
multi-chance fission; this package does not implement them, and [`WAHL_VALIDITY`](@ref) excludes
that range.

```@autodocs
Modules = [FissionFragmentsDomain]
Pages = ["charge.jl", "zp.jl", "wahl1988.jl", "wahl.jl", "domain.jl"]
```
