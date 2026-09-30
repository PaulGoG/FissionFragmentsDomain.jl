"""
    FragmentQuantity

A model quantity tabulated over fragment configurations, `q(A, Z, TKE)`, stored sparsely.
Every observable is a contraction of one of these against a weight, the isobaric charge
distribution `p(Z,A)` when reducing over charge alone and the yield `Y(A,Z,TKE)` otherwise.
Values are `Float64` or `Measurement{Float64}`, the element type stating whether the quantity
carries an uncertainty; quantities downstream of the `TXE` partition are `Float64`.
"""
struct FragmentQuantity{T <: Union{Float64, Measurement{Float64}}}
    values::Dict{Tuple{Int, Int, Float64}, T}
end

FragmentQuantity{T}() where {T} = FragmentQuantity(Dict{Tuple{Int, Int, Float64}, T}())

Base.length(q::FragmentQuantity) = length(q.values)
Base.getindex(q::FragmentQuantity, key) = q.values[key]
Base.get(q::FragmentQuantity, key, default) = get(q.values, key, default)
Base.haskey(q::FragmentQuantity, key) = haskey(q.values, key)
Base.setindex!(q::FragmentQuantity, value, key) = (q.values[key] = value)
Base.keys(q::FragmentQuantity) = keys(q.values)

function Base.show(io::IO, q::FragmentQuantity{T}) where {T}
    print(io, "FragmentQuantity{", T, "}(", length(q.values), " configurations)")
    return nothing
end

"""
    average(quantity, weights; over) -> Vector{Tuple{Float64,V}}

Contract `q(A,Z,TKE)` against a weight over the indices named in `over`, leaving the rest.

```
⟨q⟩(x) = Σ q(x,…) W(x,…) / Σ W(x,…)
```

`over` names exactly two of `(:mass, :charge, :kinetic_energy)`; the result is indexed by the
remaining one, ascending. Reducing over all three is [`total_average`](@ref). A configuration
enters both sums only where both the quantity and the weight exist, so an unsolved
configuration is omitted rather than averaged in as zero.

Uncertainties propagate by the standard expression for a weighted mean, with both terms,

```
δ²⟨q⟩ = Σᵢ [ (Wᵢ/ΣW) δqᵢ ]² + Σᵢ [ ((qᵢ − ⟨q⟩)/ΣW) δWᵢ ]²
```

the first from the quantity, the second from the weights. Both follow from `Measurement`
arithmetic, which also retains correlations between terms sharing a source, such as the
`Y(A,Z,TKE)` cells at one `(A,TKE)`, which are `p(Z,A)` times a single measurement.
"""
function average(
    quantity::FragmentQuantity{Q},
    weights::FragmentQuantity{W};
    over::NTuple{N, Symbol},
) where {Q, W, N}
    _validate_indices(over)
    length(over) == 2 || throw(
        ArgumentError(
            "average reduces to a one-dimensional distribution, so it takes exactly two " *
            "indices to reduce over; got $over",
        ),
    )
    kept = only(_kept_indices(over))
    P = typeof(one(Q) * one(W))
    numerator = Dict{Float64, P}()
    denominator = Dict{Float64, W}()

    for (key, quantity_value) in quantity.values
        weight = get(weights, key, nothing)
        weight === nothing && continue
        label = Float64(key[kept])
        numerator[label] = get(numerator, label, zero(P)) + quantity_value * weight
        denominator[label] = get(denominator, label, zero(W)) + weight
    end

    out = Tuple{Float64, P}[]
    for (label, total) in denominator
        _positive(total) || continue
        push!(out, (label, numerator[label] / total))
    end
    sort!(out; by = first)
    return out
end

"""
    total_average(quantity, weights; masses = nothing) -> Union{Float64,Measurement{Float64}}

Contract over all three indices, optionally restricting to a set of fragment masses.

`masses` separates the fragment groups: the light range gives `⟨q⟩_L`, the heavy range
`⟨q⟩_H`, and `nothing` the whole distribution. Both ranges include `A₀/2` when the domain
reaches symmetric fission, so configurations at `A₀/2` enter both groups and the two group
weights sum to slightly more than the whole. Returns
`nothing` when no configuration carries both a value and a weight.
"""
function total_average(
    quantity::FragmentQuantity{Q},
    weights::FragmentQuantity{W};
    masses = nothing,
) where {Q, W}
    P = typeof(one(Q) * one(W))
    numerator = zero(P)
    denominator = zero(W)
    counted = 0
    for (key, quantity_value) in quantity.values
        masses === nothing || key[1] in masses || continue
        weight = get(weights, key, nothing)
        weight === nothing && continue
        numerator += quantity_value * weight
        denominator += weight
        counted += 1
    end
    counted == 0 && return nothing
    _positive(denominator) || return nothing
    return numerator / denominator
end

_positive(x::Real) = x > 0
_positive(x::Measurement) = value(x) > 0

const _INDEX_NAMES = (:mass, :charge, :kinetic_energy)

function _validate_indices(over::NTuple{N, Symbol}) where {N}
    for name in over
        name in _INDEX_NAMES ||
            throw(ArgumentError("unknown index $name; expected any of $(_INDEX_NAMES)"))
    end
    allunique(over) ||
        throw(ArgumentError("indices to reduce over must be distinct: $over"))
    length(over) < 3 || throw(
        ArgumentError("reducing over every index gives a number; use `total_average`"),
    )
    return nothing
end

function _kept_indices(over::NTuple{N, Symbol}) where {N}
    kept = Int[]
    for (position, name) in pairs(_INDEX_NAMES)
        name in over || push!(kept, position)
    end
    return kept
end

"""
    weights_from(yields) -> FragmentQuantity

Use a fragment yield distribution as the weight of a contraction.
"""
weights_from(yields::FragmentYield) = FragmentQuantity(copy(yields.values))

"""
    weights_from(domain, distribution, kinetic_energies) -> FragmentQuantity

Use the isobaric charge distribution as the weight, for a run with no experimental yield:
`p(Z,A)` broadcast over `kinetic_energies`. Reducing over charge gives the same result as with
`Y = p·Y(A,TKE)`, since `Y(A,TKE)` does not depend on `Z` and cancels. Reducing over kinetic
energy does not: a flat weight centres the contraction on the midpoint of the `TKE` grid, and
every quantity built on `TXE` shifts accordingly. Where no joint `Y(A,TKE)` was measured,
[`factorized_yield`](@ref) reconstructs one from the measured marginals.
"""
function weights_from(
    domain::FragmentationDomain,
    distribution::ChargeDistribution,
    kinetic_energies,
)
    probabilities = fragment_charge_probabilities(domain, distribution)
    weights = Dict{Tuple{Int, Int, Float64}, Float64}()
    for (nuclide, p) in probabilities, TKE in kinetic_energies
        weights[(nuclide.A, nuclide.Z, Float64(TKE))] = p
    end
    return FragmentQuantity(weights)
end
