```@meta
CurrentModule = FissionFragmentsDomain
```

# Level density and energetics

The fragment level density parameter, the energy balance of a fragmentation, and the relation
between the temperature ratio of a pair and the division of its excitation energy. A
temperature-ratio extraction inverts this relation and a TXE partition applies it, so both take it
from this package.

A partition applies `E*_H/TXE = 1/(1 + R_a R_T²)` to every fragmentation with its own
`R_a = a_L/a_H`. An extraction sees only the mass-resolved multiplicity ratio `r_ν(A_H)`.
[`ChargeResolved`](@ref) inverts the partition relation exactly from `r_ν(A_H)`: the relation
reduced over the charge distribution, each fragmentation weighted by `p(Z, A_H)` and, with
[`mean_total_excitation`](@ref), by the excitation it shares out. It closes the round trip to
rounding. [`RatioOfMeans`](@ref) and [`MeanOfRatios`](@ref) replace the reduction by one effective
`R_a`; they reproduce published extractions and miss the round trip by up to 10⁻² in `R_T` at
`A_H ≈ 130`. A run records which it used in the `[domain]` table of its manifest,
[`ManifestDomain`](@ref).

```@autodocs
Modules = [FissionFragmentsDomain]
Pages = ["leveldensity.jl", "energetics.jl", "temperature_ratio.jl"]
```
