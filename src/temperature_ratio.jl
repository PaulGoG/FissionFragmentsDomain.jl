"""
    RatioAveraging

How the level density parameter ratio of complementary fragments enters the relation between
`R_T` and `E*_H/TXE` over the isobaric charge distribution at fixed `A_H`:

- [`ChargeResolved`](@ref): no effective ratio; the relation is written for every fragmentation
  and reduced over charge afterwards;
- [`RatioOfMeans`](@ref): `⟨a_L⟩/⟨a_H⟩`;
- [`MeanOfRatios`](@ref): `⟨a_L/a_H⟩`.

The last two are effective ratios. They differ from each other, and from the charge-resolved
relation, at second order in the spread of `a_L/a_H` over the charge window, and all three
coincide for a single charge. A prompt emission code that partitions each fragment pair with its
own parameters inverts exactly under `ChargeResolved` alone. Under either effective ratio the round
trip `r_ν → R_T → r_ν` misses by up to about 10⁻² in `R_T`, at the doubly magic heavy fragment.
"""
abstract type RatioAveraging end

"""
    ChargeResolved()

The relation `E*_H/TXE = 1/(1 + ρ_Z R_T²)`, with `ρ_Z = a_L/a_H` of the fragmentation
`(A_H, Z_H)`, written for every fragmentation and reduced over the charge distribution afterwards:

```
r_ν(A_H) = Σ_Z p(Z, A_H) / (1 + ρ_Z R_T²)  /  Σ_Z p(Z, A_H).
```

This is how the Point-by-Point and sequential emission treatments apply a temperature ratio: to
each fragment pair with its own level density parameters (A. Tudora, *Eur. Phys. J. A*
**58**, 126 (2022), eq. (5), doi:10.1140/epja/s10050-022-00766-y). The right-hand side falls
strictly from one to zero as `R_T` grows, so every `r_ν ∈ (0, 1)` has one root.

At the symmetric split, with a charge set invariant under `Z → Z₀ − Z`, the terms pair as `ρ` and
`1/ρ` with equal weight, and `1/(1 + ρ) + 1/(1 + 1/ρ) = 1`. So `r_ν = 1/2` returns `R_T = 1`
exactly.

    ChargeResolved(excitation)

The same relation with every fragmentation also weighted by its mean total excitation energy,
`excitation[heavy] = ⟨TXE⟩(A_H, Z_H)` (see [`mean_total_excitation`](@ref)):

```
r_ν(A_H) = Σ_Z p(Z, A_H) ⟨TXE⟩_Z / (1 + ρ_Z R_T²)  /  Σ_Z p(Z, A_H) ⟨TXE⟩_Z.
```

This is the premise of the extraction, `ν_L/ν_H = E*_L/E*_H` (Tudora and Gogita, *Eur. Phys. J.
A* **60**, 190 (2024), eq. (1), doi:10.1140/epja/s10050-024-01375-7), carried to mass-resolved
multiplicities. Each is the yield-weighted mean over charge and TKE, so a fragmentation counts in
proportion to the excitation it shares out. The Q-value and `a_L/a_H` both change across the
charge window at the `Z = 50` shell. Dropping the weight therefore moves `R_T` by up to 10⁻² at
`A_H ≈ 130`.

The two fragmentations `Z` and `Z₀ − Z` at the symmetric split are the same pair, with the same
excitation, so the identity there survives. A fragmentation absent from `excitation` is left out,
and a non-positive `⟨TXE⟩` weighs nothing: that fragmentation has no excitation to share.
"""
struct ChargeResolved <: RatioAveraging
    excitation::Union{Nothing, Dict{Nuclide, Float64}}
end

ChargeResolved() = ChargeResolved(nothing)

"""
    mean_total_excitation(masses, domain, mean_kinetic_energy) -> Dict{Nuclide,Float64}

The mean total excitation energy of every fragmentation of `domain`, keyed by its heavy
fragment, in MeV:

```
⟨TXE⟩(A_H, Z_H) = Q(A_H, Z_H) + E*_CN − ⟨TKE⟩(A_H),
```

with `Q` and the compound-nucleus excitation `E*_CN` from `masses`, and `mean_kinetic_energy`
mapping a heavy mass number to its pre-neutron `⟨TKE⟩`. `⟨TXE⟩` is linear in `TKE`, so the mean
alone carries the reduction over the kinetic-energy distribution. A heavy mass absent from
`mean_kinetic_energy`, or a fragmentation whose Q-value `masses` cannot form, is absent from the
result.
"""
function mean_total_excitation(
    masses::MassExcessTable,
    domain::FragmentationDomain,
    mean_kinetic_energy::AbstractDict{<:Integer, <:Real},
)
    E_compound = compound_nucleus_excitation(masses, domain.system)
    E_compound === nothing && throw(
        ArgumentError(
            "the compound-nucleus excitation of $(domain.system.compound) cannot be formed \
             from the mass table",
        ),
    )
    excitation = Dict{Nuclide, Float64}()
    for entry in domain.entries
        haskey(mean_kinetic_energy, entry.heavy.A) || continue
        Q = q_value(masses, domain.system, entry.heavy)
        Q === nothing && continue
        excitation[entry.heavy] =
            value(Q) + value(E_compound) - mean_kinetic_energy[entry.heavy.A]
    end
    return excitation
end

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
function level_density_ratio(::ChargeResolved, ::LevelDensityModel, ::FragmentationDomain)
    throw(
        ArgumentError(
            "ChargeResolved forms no effective level density ratio; take the relation over the \
             charge distribution with heavy_excitation_fraction or temperature_ratio at A_H",
        ),
    )
end

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

# The per-charge ratios ρ = a_L/a_H at A_H and their normalised weights, p(Z, A_H) or
# p(Z, A_H)⟨TXE⟩, over the fragmentations whose parameters and excitation exist; `nothing`
# where none does.
function _charge_terms(
    averaging::ChargeResolved,
    model::LevelDensityModel,
    domain::FragmentationDomain,
    A_H::Integer,
)
    ρ = Float64[]
    w = Float64[]
    excitation = averaging.excitation
    for entry in domain.entries
        entry.heavy.A == A_H || continue
        a_H = level_density_parameter(model, entry.heavy)
        a_L = level_density_parameter(model, entry.light)
        (a_H === nothing || a_L === nothing) && continue
        weight = entry.probability
        if excitation !== nothing
            haskey(excitation, entry.heavy) || continue
            weight *= max(excitation[entry.heavy], 0.0)
        end
        push!(ρ, a_L / a_H)
        push!(w, weight)
    end
    total = sum(w; init = 0.0)
    total > 0 || return nothing
    return ρ, w ./ total
end

_resolved_fraction(ρ, w, R_T) = sum(w[i] / (1 + ρ[i] * R_T^2) for i in eachindex(ρ))

"""
    heavy_excitation_fraction(averaging, model, domain, A_H, R_T) -> Union{Float64,Nothing}

`E*_H/TXE` at the heavy mass `A_H` of `domain`, reduced over its charge distribution as
`averaging` prescribes: the charge-weighted mean of the per-fragmentation relation for
[`ChargeResolved`](@ref), the relation at the effective ratio for [`RatioOfMeans`](@ref) and
[`MeanOfRatios`](@ref). Level density parameters come from `model`. Returns `nothing` where no
fragmentation at `A_H` has both parameters.
"""
function heavy_excitation_fraction(
    averaging::ChargeResolved,
    model::LevelDensityModel,
    domain::FragmentationDomain,
    A_H::Integer,
    R_T::Real,
)
    R_T > 0 || throw(DomainError(R_T, "R_T must be positive, got $R_T"))
    terms = _charge_terms(averaging, model, domain, A_H)
    terms === nothing && return nothing
    return _resolved_fraction(terms..., R_T)
end

function heavy_excitation_fraction(
    averaging::RatioAveraging,
    model::LevelDensityModel,
    domain::FragmentationDomain,
    A_H::Integer,
    R_T::Real,
)
    R_a = _charge_averaged_ratio(averaging, model, domain, A_H)
    R_a === nothing && return nothing
    return heavy_excitation_fraction(R_T, R_a)
end

"""
    temperature_ratio(averaging, model, domain, A_H, r_ν) -> Union{Float64,Nothing}

The temperature ratio at the heavy mass `A_H` that reproduces the multiplicity ratio `r_ν`, the
inverse of [`heavy_excitation_fraction`](@ref) with the same arguments. For
[`ChargeResolved`](@ref) the root is taken by bisection. It is bracketed by the closed-form roots
at the largest and the smallest `a_L/a_H` of the charge window, and resolved to rounding. For an
effective ratio it is the closed form. Returns `nothing` where no fragmentation at `A_H` has both
parameters.

# Examples

```jldoctest
julia> table = read_mass_excess_table(String(AME2020_MASS_EXCESS_FILE));

julia> system = spontaneous_fission(Nuclide(98, 252));

julia> domain = fragmentation_domain(system, WAHL_1988[(98, 252)], 126:174);

julia> temperature_ratio(ChargeResolved(), BackShiftedFermiGas(table), domain, 126, 0.5) ≈ 1
true
```
"""
function temperature_ratio(
    averaging::ChargeResolved,
    model::LevelDensityModel,
    domain::FragmentationDomain,
    A_H::Integer,
    r_ν::Real,
)
    0 < r_ν < 1 || throw(DomainError(r_ν, "r_ν must lie in (0, 1), got $r_ν"))
    terms = _charge_terms(averaging, model, domain, A_H)
    terms === nothing && return nothing
    ρ, w = terms
    # The fraction lies between the single-charge relations at the extreme ratios, so their
    # closed-form roots bracket the root.
    lower = temperature_ratio(r_ν, maximum(ρ))
    upper = temperature_ratio(r_ν, minimum(ρ))
    lower == upper && return lower
    for _ in 1:200
        middle = (lower + upper) / 2
        (middle == lower || middle == upper) && break
        if _resolved_fraction(ρ, w, middle) > r_ν
            lower = middle
        else
            upper = middle
        end
    end
    # Of the two adjacent floating-point bounds, the one closer to the target.
    return abs(_resolved_fraction(ρ, w, lower) - r_ν) <=
           abs(_resolved_fraction(ρ, w, upper) - r_ν) ? lower : upper
end

function temperature_ratio(
    averaging::RatioAveraging,
    model::LevelDensityModel,
    domain::FragmentationDomain,
    A_H::Integer,
    r_ν::Real,
)
    R_a = _charge_averaged_ratio(averaging, model, domain, A_H)
    R_a === nothing && return nothing
    return temperature_ratio(r_ν, R_a)
end

"""
    temperature_ratio_slope(averaging, model, domain, A_H, R_T) -> Union{Float64,Nothing}

`∂R_T/∂r_ν` at the heavy mass `A_H` and temperature ratio `R_T`: the factor by which
[`temperature_ratio`](@ref) with the same arguments carries an uncertainty or a covariance of
`r_ν` into `R_T`. For [`ChargeResolved`](@ref) it is the reciprocal of

```
∂r_ν/∂R_T = −Σ_Z p(Z, A_H) 2 ρ_Z R_T / (1 + ρ_Z R_T²)²  /  Σ_Z p(Z, A_H).
```

Returns `nothing` where no fragmentation at `A_H` has both parameters.
"""
function temperature_ratio_slope(
    averaging::ChargeResolved,
    model::LevelDensityModel,
    domain::FragmentationDomain,
    A_H::Integer,
    R_T::Real,
)
    R_T > 0 || throw(DomainError(R_T, "R_T must be positive, got $R_T"))
    terms = _charge_terms(averaging, model, domain, A_H)
    terms === nothing && return nothing
    ρ, w = terms
    derivative = -sum(w[i] * 2 * ρ[i] * R_T / (1 + ρ[i] * R_T^2)^2 for i in eachindex(ρ))
    return 1 / derivative
end

function temperature_ratio_slope(
    averaging::RatioAveraging,
    model::LevelDensityModel,
    domain::FragmentationDomain,
    A_H::Integer,
    R_T::Real,
)
    R_a = _charge_averaged_ratio(averaging, model, domain, A_H)
    R_a === nothing && return nothing
    return temperature_ratio_slope(R_T, R_a, heavy_excitation_fraction(R_T, R_a))
end

"""
    RATIO_AVERAGINGS

The spellings of the ratio averagings: `"charge_resolved"` for [`ChargeResolved`](@ref),
`"ratio_of_means"` for [`RatioOfMeans`](@ref) and `"mean_of_ratios"` for
[`MeanOfRatios`](@ref).
"""
const RATIO_AVERAGINGS = ("charge_resolved", "ratio_of_means", "mean_of_ratios")

"""
    ratio_averaging_label(averaging::RatioAveraging) -> String
    ratio_averaging(label::AbstractString) -> RatioAveraging

The spelling of a [`RatioAveraging`](@ref) in configurations and run records, one of
[`RATIO_AVERAGINGS`](@ref), and its inverse. An unknown label is an `ArgumentError`.
"""
ratio_averaging_label(::ChargeResolved) = "charge_resolved"
ratio_averaging_label(::RatioOfMeans) = "ratio_of_means"
ratio_averaging_label(::MeanOfRatios) = "mean_of_ratios"

@doc (@doc ratio_averaging_label)
function ratio_averaging(label::AbstractString)
    label == "charge_resolved" && return ChargeResolved()
    label == "ratio_of_means" && return RatioOfMeans()
    label == "mean_of_ratios" && return MeanOfRatios()
    throw(ArgumentError("unknown ratio averaging $(repr(label)); one of \
             $(join(map(repr, RATIO_AVERAGINGS), ", "))"))
end

"""
    level_density_label(model::LevelDensityModel) -> String

The spelling of a level density model in configurations and run records, one of
[`LEVEL_DENSITY_MODELS`](@ref).
"""
level_density_label(::BackShiftedFermiGas) = "BSFG"
level_density_label(::GilbertCameron) = "GC"
