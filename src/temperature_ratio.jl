"""
    RatioAveraging

The order in which the level density parameter ratio of complementary fragments is reduced over
the isobaric charge distribution at fixed `A_H`: [`RatioOfMeans`](@ref), `⟨a_L⟩/⟨a_H⟩`, or
[`MeanOfRatios`](@ref), `⟨a_L/a_H⟩`. The two differ by Jensen's inequality and coincide only for
a degenerate charge distribution. An extraction of `R_T(A_H)` and the partition that applies it
must use the same order, or the round trip `ν(A) → R_T → ν(A)` does not close.
"""
abstract type RatioAveraging end

"""
    RatioOfMeans()

`R_a = ⟨a_L⟩/⟨a_H⟩`: the Fermi-gas relation `⟨E*⟩ = ⟨a⟩ T²` for the mass-resolved fragment, with a
temperature that depends on mass alone. At the symmetric split, with a charge set invariant under
`Z → Z₀ − Z`, both averages run over the same nuclides and `R_a = 1` exactly.
"""
struct RatioOfMeans <: RatioAveraging end

"""
    MeanOfRatios()

`R_a = ⟨a_L/a_H⟩`: the partition written for each fragmentation `(A_H, Z_H)` and averaged
afterwards. At the symmetric split the per-charge ratios come in reciprocal pairs of equal
weight, so this order gives `R_a > 1` there by a part in a thousand.
"""
struct MeanOfRatios <: RatioAveraging end

"""
    level_density_ratio(averaging, model, domain) -> Dict{Int,Float64}

The level density parameter ratio of complementary fragments, `R_a(A_H) = a_L/a_H`, averaged over
the isobaric charge distribution in the order `averaging` fixes, at every heavy mass of the
domain. Each fragmentation is weighted by its `p(Z,A)`; a fragmentation where either parameter is
undefined is left out of both averages, and the weights are renormalised over those that remain.
A heavy mass with no defined pair, or a non-positive result, is absent from the result.

At `A_H = A₀/2` the heavy entries run over the retained window; when it is its own mirror under
`Z → Z₀ − Z` (see [`symmetric_charge_set_is_invariant`](@ref)) every split enters in both
labellings with equal weight.
"""
function level_density_ratio(
    averaging::RatioAveraging,
    model::LevelDensityModel,
    domain::FragmentationDomain,
)
    ratio = Dict{Int, Float64}()
    for A_H in domain.heavy_masses
        R_a = _charge_averaged_ratio(averaging, model, domain, A_H)
        R_a === nothing || R_a <= 0 || (ratio[A_H] = R_a)
    end
    return ratio
end

function _charge_averaged_ratio(
    ::MeanOfRatios,
    model::LevelDensityModel,
    domain::FragmentationDomain,
    A_H::Integer,
)
    numerator = 0.0
    weights = 0.0
    for entry in domain.entries
        entry.heavy.A == A_H || continue
        a_H = level_density_parameter(model, entry.heavy)
        a_L = level_density_parameter(model, entry.light)
        (a_H === nothing || a_L === nothing) && continue
        numerator += entry.probability * a_L / a_H
        weights += entry.probability
    end
    return weights > 0 ? numerator / weights : nothing
end

function _charge_averaged_ratio(
    ::RatioOfMeans,
    model::LevelDensityModel,
    domain::FragmentationDomain,
    A_H::Integer,
)
    light = 0.0
    heavy = 0.0
    for entry in domain.entries
        entry.heavy.A == A_H || continue
        a_H = level_density_parameter(model, entry.heavy)
        a_L = level_density_parameter(model, entry.light)
        (a_H === nothing || a_L === nothing) && continue
        light += entry.probability * a_L
        heavy += entry.probability * a_H
    end
    return heavy > 0 ? light / heavy : nothing
end

"""
    heavy_excitation_fraction(R_T, R_a) -> Float64

The share of the total excitation energy the heavy fragment carries, given the temperature ratio
`R_T = T_L/T_H` and the level density parameter ratio `R_a = a_L/a_H` of the pair:

```
E*_H/TXE = 1 / (1 + R_a R_T²),
```

from `E* = a T²` for each fragment, so that `E*_L/E*_H = R_a R_T²`. The inverse, from a
multiplicity ratio, is [`temperature_ratio`](@ref).
"""
function heavy_excitation_fraction(R_T::Real, R_a::Real)
    (R_T > 0 && R_a > 0) || throw(
        DomainError((R_T, R_a), "R_T and R_a must be positive, got R_T = $R_T, R_a = $R_a"),
    )
    return 1 / (1 + R_a * R_T^2)
end

"""
    temperature_ratio(r_ν, R_a) -> Float64

The temperature ratio of complementary fully accelerated fragments,

```
R_T = T_L/T_H = [(1 − r_ν) / (R_a r_ν)]^(1/2),
```

from the multiplicity ratio `r_ν = ν_H/(ν_L + ν_H)`, identified with `E*_H/TXE`, and the level
density parameter ratio `R_a = a_L/a_H`. It inverts [`heavy_excitation_fraction`](@ref) exactly
when both use the same `R_a`; see [`RatioAveraging`](@ref).

# Examples

```jldoctest
julia> round(temperature_ratio(0.40, 1.05); digits = 4)
1.1952

julia> heavy_excitation_fraction(temperature_ratio(0.40, 1.05), 1.05) ≈ 0.40
true
```
"""
function temperature_ratio(r_ν::Real, R_a::Real)
    (0 < r_ν < 1 && R_a > 0) || throw(
        DomainError(
            (r_ν, R_a),
            "r_ν must lie in (0, 1) and R_a be positive, got r_ν = $r_ν, R_a = $R_a",
        ),
    )
    return sqrt((1 - r_ν) / (R_a * r_ν))
end

"""
    temperature_ratio_slope(R_T, R_a, r_ν) -> Float64

`∂R_T/∂r_ν = −1 / (2 R_T R_a r_ν²)`, the factor by which [`temperature_ratio`](@ref) carries an
uncertainty or a covariance of `r_ν` into `R_T`; `R_a` carries none of its own.
"""
temperature_ratio_slope(R_T::Real, R_a::Real, r_ν::Real) = -1 / (2 * R_T * R_a * r_ν^2)
