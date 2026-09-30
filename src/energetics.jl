"""
    compound_nucleus_excitation(table, system) -> Union{Measurement{Float64},Nothing}

Excitation energy of the compound nucleus in MeV.

Zero for spontaneous fission. For neutron-induced fission it is the neutron separation energy
of the compound nucleus plus the incident kinetic energy, `E*_CN = Sₙ(A₀,Z₀) + Eₙ`. Returns
`nothing` when `Sₙ` cannot be formed from the mass table.
"""
function compound_nucleus_excitation(table::MassExcessTable, system::FissioningSystem)
    is_spontaneous(system) && return measurement(0.0, 0.0)
    Sₙ = neutron_separation_energy(table, system.compound)
    Sₙ === nothing && return nothing
    return Sₙ + system.incident_energy
end

"""
    total_excitation_energy(Q, E_compound, TKE) -> Union{Measurement{Float64},Nothing}

Total excitation energy shared by the complementary fragments at full acceleration, in MeV:

```
TXE = Q + E*_CN − TKE
```

Returns `nothing` when `TXE ≤ 0`, an unphysical configuration; a deterministic sweep over TKE
therefore produces fewer records than its declared grid. `TKE` carries no uncertainty, being the
swept independent variable of the model.

Uncertainty propagation stops here. `TXE` carries the uncertainty of the mass excesses, but the
TXE partition methods that consume it — modelling at scission, `E*_H/TXE (A_H)`, `R_T (A_H)` —
are fitted curves without quoted uncertainty, so `E*_L` and `E*_H` inherit none; a propagated
`σ` would omit the dominant model term.
"""
function total_excitation_energy(
    Q::Measurement{Float64},
    E_compound::Measurement{Float64},
    TKE::Real,
)
    TXE = Q + E_compound - TKE
    return value(TXE) > 0 ? TXE : nothing
end

"""
    kinetic_energy(system, fragment, TKE) -> Float64

Kinetic energy of one fragment in MeV, from conservation of momentum in a binary split:

```
KE(A) = TKE (A₀ − A) / A₀
```
"""
function kinetic_energy(system::FissioningSystem, fragment::Nuclide, TKE::Real)
    return TKE * (system.compound.A - fragment.A) / system.compound.A
end

"""
    Energetics

The energy balance of one fragmentation at one `TKE`, carrying the quantities that precede the
TXE partition.
"""
struct Energetics
    fragmentation::Fragmentation
    TKE::Float64
    Q::Measurement{Float64}
    TXE::Measurement{Float64}

    function Energetics(
        fragmentation::Fragmentation,
        TKE::Real,
        Q::Measurement{Float64},
        TXE::Measurement{Float64},
    )
        return new(fragmentation, Float64(TKE), Q, TXE)
    end
end

function Base.show(io::IO, energetics::Energetics)
    print(
        io,
        "Energetics(A_H=",
        energetics.fragmentation.heavy.A,
        ", Z_H=",
        energetics.fragmentation.heavy.Z,
        ", TKE=",
        energetics.TKE,
        ", TXE=",
        energetics.TXE,
        ")",
    )
    return nothing
end

"""
    energetics(table, system, fragmentation, TKE) -> Union{Energetics,Nothing}

The energy balance of one fragmentation at one kinetic energy, or `nothing` when the
configuration is unphysical — a missing mass excess, a non-positive `Q`, or `TXE ≤ 0`.
"""
function energetics(
    table::MassExcessTable,
    system::FissioningSystem,
    fragmentation::Fragmentation,
    TKE::Real,
)
    Q = q_value(table, system, fragmentation.heavy)
    Q === nothing && return nothing
    E_compound = compound_nucleus_excitation(table, system)
    E_compound === nothing && return nothing
    TXE = total_excitation_energy(Q, E_compound, TKE)
    TXE === nothing && return nothing
    return Energetics(fragmentation, Float64(TKE), Q, TXE)
end

"""
    Exclusion

Why one `(fragmentation, TKE)` configuration was left out of a sweep, recorded so that the
coverage of the declared grid can be audited.
"""
struct Exclusion
    fragmentation::Fragmentation
    TKE::Float64
    reason::Symbol

    function Exclusion(fragmentation::Fragmentation, TKE::Real, reason::Symbol)
        return new(fragmentation, Float64(TKE), reason)
    end
end

"""
    _validate_kinetic_energies(kinetic_energies)

Check the swept `TKE` grid as [`fragmentation_domain`](@ref) checks the mass range: non-empty,
finite, positive and strictly ascending. A negative `TKE` raises `TXE = Q + E*_CN − TKE` and
would pass the `TXE > 0` test; an empty grid would give a sweep of zero records,
indistinguishable from physical exclusion; and the sweep emits in grid order.
"""
function _validate_kinetic_energies(kinetic_energies::AbstractVector{<:Real})
    isempty(kinetic_energies) &&
        throw(ArgumentError("the kinetic energy grid is empty; nothing would be swept"))
    # One adjacency test covers ascent and distinctness and avoids Base's `allunique` path for
    # abstractly typed vectors, which static analysis cannot follow through reinterpret/reshape.
    start = firstindex(kinetic_energies)
    for index in eachindex(kinetic_energies)
        TKE = Float64(kinetic_energies[index])
        isfinite(TKE) ||
            throw(ArgumentError("kinetic energy number $index of the grid is not finite"))
        TKE > 0 || throw(
            ArgumentError(
                "total kinetic energy must be positive; entry $index is $TKE MeV",
            ),
        )
        index == start && continue
        previous = Float64(kinetic_energies[index - 1])
        previous < TKE || throw(
            ArgumentError(
                "the grid must ascend strictly; entry $index is $TKE after $previous",
            ),
        )
    end
    return nothing
end

"""
    sweep_energetics(table, system, domain, kinetic_energies)
        -> (Vector{Energetics}, Vector{Exclusion})

Evaluate the energy balance over the whole domain and TKE grid, returning both what survived
and what did not, each exclusion carrying its reason: `:missing_mass`, `:nonpositive_q` or
`:nonpositive_txe`.
"""
function sweep_energetics(
    table::MassExcessTable,
    system::FissioningSystem,
    domain::FragmentationDomain,
    kinetic_energies::AbstractVector{<:Real},
)
    _validate_kinetic_energies(kinetic_energies)

    accepted = Energetics[]
    excluded = Exclusion[]
    sizehint!(accepted, length(domain) * length(kinetic_energies))

    E_compound = compound_nucleus_excitation(table, system)
    E_compound === nothing && throw(
        ArgumentError(
            "the neutron separation energy of the compound nucleus $(system.compound) is not " *
            "available from $(table.source); the excitation energy cannot be formed",
        ),
    )

    for fragmentation in domain
        Q = q_value(table, system, fragmentation.heavy)
        if Q === nothing
            reason = _q_failure_reason(table, system, fragmentation.heavy)
            for TKE in kinetic_energies
                push!(excluded, Exclusion(fragmentation, Float64(TKE), reason))
            end
            continue
        end
        for TKE in kinetic_energies
            TXE = total_excitation_energy(Q, E_compound, TKE)
            if TXE === nothing
                push!(excluded, Exclusion(fragmentation, Float64(TKE), :nonpositive_txe))
            else
                push!(accepted, Energetics(fragmentation, Float64(TKE), Q, TXE))
            end
        end
    end
    return accepted, excluded
end

function _q_failure_reason(table::MassExcessTable, system::FissioningSystem, heavy::Nuclide)
    light = complementary_fragment(system, heavy)
    if !haskey(table, heavy) || !haskey(table, light) || !haskey(table, system.compound)
        return :missing_mass
    end
    return :nonpositive_q
end

"""
    exclusion_summary(excluded) -> Dict{Symbol,Int}

Count the excluded configurations by reason, for the run's provenance record.
"""
function exclusion_summary(excluded::AbstractVector{Exclusion})
    counts = Dict{Symbol, Int}()
    for exclusion in excluded
        counts[exclusion.reason] = get(counts, exclusion.reason, 0) + 1
    end
    return counts
end

"""
    fragment_energy_per_nucleon(fragmentation, TKE) -> (light, heavy)

The kinetic energy per nucleon of each fully accelerated fragment, in MeV.

Equal and opposite momenta give `KE_L = TKE·A_H/A₀` and `KE_H = TKE·A_L/A₀`; dividing each by
its own mass,

```
E_f^L = (A_H/A_L) TKE/A₀,    E_f^H = (A_L/A_H) TKE/A₀
```

so `E_f^L > E_f^H` and the light-fragment spectrum is the broader. The neutron mass is taken as
one atomic mass unit; `m_n/m_u = 1.00867` would raise both by 0.9 %. `E_f` is constant along an
emission chain: under the recoil-free scaling `KE_k = KE (A−k)/A` the ratio `KE_k/(A−k)` is
invariant.
"""
function fragment_energy_per_nucleon(fragmentation::Fragmentation, TKE::Real)
    A₀ = fragmentation.light.A + fragmentation.heavy.A
    light = (fragmentation.heavy.A / fragmentation.light.A) * TKE / A₀
    heavy = (fragmentation.light.A / fragmentation.heavy.A) * TKE / A₀
    return (light, heavy)
end
