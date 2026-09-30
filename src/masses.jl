"""
    MassExcessTable

Tabulated mass excesses with their uncertainties, in MeV, indexed by nuclide. Values are
`Measurement`s, so derived quantities carry their uncertainty and a nuclide entering an
expression twice is correlated with itself.
"""
struct MassExcessTable
    Δ::Dict{Nuclide, Measurement{Float64}}
    source::String

    function MassExcessTable(Δ::Dict{Nuclide, Measurement{Float64}}, source::AbstractString)
        return new(Δ, String(source))
    end
end

Base.length(table::MassExcessTable) = length(table.Δ)
Base.haskey(table::MassExcessTable, nuclide::Nuclide) = haskey(table.Δ, nuclide)

function Base.show(io::IO, table::MassExcessTable)
    print(
        io,
        "MassExcessTable(",
        length(table.Δ),
        " nuclides from ",
        basename(table.source),
        ")",
    )
    return nothing
end

"""
    MASS_EXCESS_SPEC

Layout of the tabulated mass-excess files: one header row, then
`Z A element mass_excess mass_excess_uncertainty`, whitespace-delimited, with the mass excess and
its uncertainty in keV. This is the layout of the shipped file
`data/reference/mass_excess_ame2020.dat`.
"""
const MASS_EXCESS_SPEC =
    TableSpec([:Z, :A, :symbol, :mass_excess, :mass_excess_uncertainty]; skip = 1)

"""
    HEADERLESS_MASS_EXCESS_SPEC

The same five columns without a header row.

```julia
read_mass_excess_table(path; spec = HEADERLESS_MASS_EXCESS_SPEC)
```
"""
const HEADERLESS_MASS_EXCESS_SPEC =
    TableSpec([:Z, :A, :symbol, :mass_excess, :mass_excess_uncertainty]; skip = 0)

const KEV_TO_MEV = 1.0e-3

"""
    read_mass_excess_table(path; spec = MASS_EXCESS_SPEC) -> MassExcessTable

Read a mass-excess table, converting from keV to MeV on load so that MeV is the only energy
unit inside the package. Columns are taken by position, not by header text: charge, mass,
element symbol, mass excess, uncertainty. Throws if the file repeats a nuclide.
"""
function read_mass_excess_table(path::AbstractString; spec::TableSpec = MASS_EXCESS_SPEC)
    table = read_delimited_table(path, spec)
    protons = integer_column(table, :Z)
    masses = integer_column(table, :A)
    excess = column(table, :mass_excess, Float64)
    uncertainty = column(table, :mass_excess_uncertainty, Float64)

    entries = Dict{Nuclide, Measurement{Float64}}()
    sizehint!(entries, length(table))
    for row in eachindex(protons)
        nuclide = Nuclide(protons[row], masses[row])
        haskey(entries, nuclide) && throw(
            ArgumentError(
                "$(table.path) row $row: nuclide $nuclide appears more than once",
            ),
        )
        entries[nuclide] =
            measurement(excess[row] * KEV_TO_MEV, uncertainty[row] * KEV_TO_MEV)
    end
    return MassExcessTable(entries, String(path))
end

"""
    mass_excess(table, nuclide) -> Union{Measurement{Float64},Nothing}

Mass excess in MeV, or `nothing` when the nuclide is absent from the table. A deterministic
fragmentation domain reaches nuclides with no measured or extrapolated mass; those
configurations are excluded rather than treated as errors.
"""
function mass_excess(table::MassExcessTable, nuclide::Nuclide)
    return get(table.Δ, nuclide, nothing)
end

"""
    separation_energy(table, nuclide, particle) -> Union{Measurement{Float64},Nothing}

Energy in MeV required to remove `particle` from `nuclide`,

```
S = Δ(A − Aₚ, Z − Zₚ) + Δ(Aₚ, Zₚ) − Δ(A, Z)
```

returning `nothing` if any of the three mass excesses is unavailable.
"""
function separation_energy(table::MassExcessTable, nuclide::Nuclide, particle::Nuclide)
    particle.A <= nuclide.A && particle.Z <= nuclide.Z ||
        throw(ArgumentError("cannot remove $particle from the smaller nuclide $nuclide"))
    residual = Nuclide(nuclide.Z - particle.Z, nuclide.A - particle.A)
    Δ_residual = mass_excess(table, residual)
    Δ_particle = mass_excess(table, particle)
    Δ_nuclide = mass_excess(table, nuclide)
    (Δ_residual === nothing || Δ_particle === nothing || Δ_nuclide === nothing) &&
        return nothing
    return Δ_residual + Δ_particle - Δ_nuclide
end

"""
    neutron_separation_energy(table, nuclide) -> Union{Measurement{Float64},Nothing}

Neutron separation energy `Sₙ` in MeV.
"""
neutron_separation_energy(table::MassExcessTable, nuclide::Nuclide) =
    separation_energy(table, nuclide, NEUTRON)

"""
    q_value(table, system, heavy) -> Union{Measurement{Float64},Nothing}

Energy released at scission for the split of the compound nucleus into `heavy` and its
complement, in MeV:

```
Q = Δ(A₀, Z₀) − [Δ(A_H, Z_H) + Δ(A_L, Z_L)]
```

Returns `nothing` when a mass excess is unavailable, and also when `Q ≤ 0`, which is not a
physical fission channel.
"""
function q_value(table::MassExcessTable, system::FissioningSystem, heavy::Nuclide)
    light = complementary_fragment(system, heavy)
    Δ_compound = mass_excess(table, system.compound)
    Δ_heavy = mass_excess(table, heavy)
    Δ_light = mass_excess(table, light)
    (Δ_compound === nothing || Δ_heavy === nothing || Δ_light === nothing) && return nothing
    Q = Δ_compound - (Δ_heavy + Δ_light)
    return value(Q) > 0 ? Q : nothing
end

"""
    complementary_fragment(system, fragment) -> Nuclide

The partner fragment of a binary split, fixed by conservation of nucleon number and charge:
`A_L = A₀ − A_H`, `Z_L = Z₀ − Z_H`.
"""
function complementary_fragment(system::FissioningSystem, fragment::Nuclide)
    fragment.A < system.compound.A && fragment.Z < system.compound.Z || throw(
        ArgumentError(
            "$fragment cannot be a fragment of $(system.compound): a binary split leaves both " *
            "partners strictly lighter and less charged than the compound nucleus",
        ),
    )
    return Nuclide(system.compound.Z - fragment.Z, system.compound.A - fragment.A)
end
