```@meta
CurrentModule = FissionFragmentsDomain
```

# Level density and energetics

The fragment level density parameter, the energy balance of a fragmentation, and the relation
between the temperature ratio of a pair and the division of its excitation energy — the relation
a temperature-ratio extraction inverts and a TXE partition applies, which is why both sides take it
from here.

A partition applies `E*_H/TXE = 1/(1 + ρ R_T²)` to every fragmentation with its own
`ρ = a_L/a_H`. An extraction sees only the mass-resolved multiplicity ratio `r_ν(A_H)`. The
exact inverse of the one through the other is [`ChargeResolved`](@ref): the relation reduced over
the charge distribution, each fragmentation weighted by `p(Z, A_H)` and, with
[`mean_total_excitation`](@ref), by the excitation it shares out. It closes the round trip to
rounding. [`RatioOfMeans`](@ref) and [`MeanOfRatios`](@ref) replace the reduction by one effective
`ρ`; they reproduce published extractions and miss the round trip by up to 10⁻² in `R_T` at
`A_H ≈ 130`. A run records which it used in the `[domain]` table of its manifest,
[`ManifestDomain`](@ref).

```@autodocs
Modules = [FissionFragmentsDomain]
Pages = ["leveldensity.jl", "energetics.jl", "temperature_ratio.jl"]
```
