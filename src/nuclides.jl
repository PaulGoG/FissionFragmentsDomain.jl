"""
    Nuclide(Z, A)

A nuclide identified by its proton number `Z` and mass number `A`.
"""
struct Nuclide
    Z::Int
    A::Int

    function Nuclide(Z::Integer, A::Integer)
        Z >= 0 || throw(ArgumentError("proton number must be non-negative, got Z = $Z"))
        A >= Z || throw(
            ArgumentError(
                "mass number must be at least the proton number, got A = $A, Z = $Z",
            ),
        )
        return new(Int(Z), Int(A))
    end
end

"""
    neutron_number(nuclide) -> Int

Neutron number `N = A - Z`.
"""
neutron_number(nuclide::Nuclide) = nuclide.A - nuclide.Z

Base.show(io::IO, nuclide::Nuclide) =
    print(io, "Nuclide(Z=", nuclide.Z, ", A=", nuclide.A, ")")

const NEUTRON = Nuclide(0, 1)

"""
    CHANNELS

Entrance channels of a fissioning system: spontaneous fission and thermal, resonance-region and
fast neutron-induced fission. The channel distinguishes systems that share the reaction code
`n,f` and names the system in a path; see [`system_label`](@ref).
"""
const CHANNELS = ("sf", "nth", "nres", "nfast")

"""
    CHANNEL_REACTION

Reaction code of each entrance channel, derived from the channel rather than configured.
"""
const CHANNEL_REACTION =
    Dict("sf" => "0,f", "nth" => "n,f", "nres" => "n,f", "nfast" => "n,f")

"""
    REACTIONS

Number of neutrons absorbed by the target in each reaction, giving the compound mass
`A₀ = A_target + REACTIONS[reaction]`.
"""
const REACTIONS = Dict("0,f" => 0, "n,f" => 1)

"""
    FissioningSystem(compound, target, incident_energy, channel)

The fissioning system, holding the target and the compound nucleus as separate fields. For
spontaneous fission the two coincide and `incident_energy` is zero; for neutron-induced fission
the compound nucleus is the target plus one neutron, and `incident_energy` is the kinetic energy
of the incident neutron in MeV. The compound nucleus defines the fragmentation domain; the target
names the fission case.

# Examples

```jldoctest
julia> spontaneous_fission(Nuclide(98, 252))
FissioningSystem(252Cf(sf))

julia> neutron_induced_fission(Nuclide(94, 239), 2.53e-8, "nth")
FissioningSystem(239Pu(nth,f) at E = 2.53e-8 MeV)
```
"""
struct FissioningSystem
    compound::Nuclide
    target::Nuclide
    incident_energy::Float64
    channel::String

    function FissioningSystem(
        compound::Nuclide,
        target::Nuclide,
        incident_energy::Real,
        channel::AbstractString,
    )
        channel in CHANNELS || throw(ArgumentError("the entrance channel must be one of \
                 $(join(map(repr, CHANNELS), ", ")), got $(repr(channel))"))
        incident_energy >= 0 || throw(
            ArgumentError("incident energy must be non-negative, got $incident_energy MeV"),
        )
        if channel == "sf"
            compound == target || throw(
                ArgumentError(
                    "spontaneous fission requires the compound and target nuclides to coincide, " *
                    "got compound $compound and target $target",
                ),
            )
            iszero(incident_energy) || throw(
                ArgumentError(
                    "spontaneous fission requires zero incident energy, got $incident_energy MeV",
                ),
            )
        end
        return new(compound, target, Float64(incident_energy), String(channel))
    end
end

"""
    is_spontaneous(system) -> Bool

Whether the system fissions spontaneously, which is the `"sf"` entrance channel and no other.
"""
is_spontaneous(system::FissioningSystem) = system.channel == "sf"

"""
    reaction(system) -> String

The reaction code of the system's entrance channel, `"0,f"` or `"n,f"`.
"""
reaction(system::FissioningSystem) = CHANNEL_REACTION[system.channel]

"""
    spontaneous_fission(nuclide) -> FissioningSystem

The spontaneously fissioning system formed by `nuclide`.
"""
spontaneous_fission(nuclide::Nuclide) = FissioningSystem(nuclide, nuclide, 0.0, "sf")

"""
    neutron_induced_fission(target, incident_energy, channel) -> FissioningSystem

The system formed by `target` capturing a neutron of kinetic energy `incident_energy` in MeV,
on the entrance channel `channel` — one of `"nth"`, `"nres"` or `"nfast"`.

The compound nucleus carries `REACTIONS["n,f"]` more mass units than the target. Consistency of
the channel with the incident energy is not checked, the thermal, resonance and fast boundaries
being conventions rather than thresholds.
"""
function neutron_induced_fission(
    target::Nuclide,
    incident_energy::Real,
    channel::AbstractString,
)
    channel == "sf" && throw(
        ArgumentError("\"sf\" is spontaneous fission; use `spontaneous_fission` instead"),
    )
    compound = Nuclide(target.Z, target.A + REACTIONS["n,f"])
    return FissioningSystem(compound, target, incident_energy, channel)
end

"""
    system_label(system) -> String

The one token that names a fissioning system, `<ElementSymbol><A><_channel>`: `Cf252_sf`,
`U235_nth`, `U235_nres`. It names the output directory of a run, its input directory, and its
configuration file. The typeset form is [`system_notation`](@ref).

# Examples

```jldoctest
julia> system_label(spontaneous_fission(Nuclide(98, 252)))
"Cf252_sf"

julia> system_label(neutron_induced_fission(Nuclide(92, 235), 5.8013e-4, "nres"))
"U235_nres"
```
"""
system_label(system::FissioningSystem) =
    string(element_symbol(system.target.Z), system.target.A, "_", system.channel)

"""
    system_notation(system) -> String

The typeset form of a fissioning system for a figure or a caption: `²⁵²Cf(sf)`, `²³⁵U(nth,f)`.

The channel appears in its plain spelling inside the reaction parentheses, following the
notation of the fission literature.

# Examples

```jldoctest
julia> system_notation(spontaneous_fission(Nuclide(98, 252)))
"²⁵²Cf(sf)"
```
"""
function system_notation(system::FissioningSystem)
    superscripts = Dict(
        '0' => '⁰',
        '1' => '¹',
        '2' => '²',
        '3' => '³',
        '4' => '⁴',
        '5' => '⁵',
        '6' => '⁶',
        '7' => '⁷',
        '8' => '⁸',
        '9' => '⁹',
    )
    mass = map(digit -> superscripts[digit], string(system.target.A))
    inside = is_spontaneous(system) ? "sf" : string(system.channel, ",f")
    return string(mass, element_symbol(system.target.Z), "(", inside, ")")
end

"""
    system_record(system) -> Dict{String,Any}

The fissioning system as a TOML table: `label`, `notation`, `target_A`, `target_Z`, `channel`,
`reaction`, `incident_energy_MeV`, `compound_A`, `compound_Z`. It is the `[system]` table of a
temperature-ratio run record and of run metadata, written by producers and checked by consumers
against their own system.
"""
function system_record(system::FissioningSystem)
    return Dict{String, Any}(
        "label" => system_label(system),
        "notation" => system_notation(system),
        "target_A" => system.target.A,
        "target_Z" => system.target.Z,
        "channel" => system.channel,
        "reaction" => reaction(system),
        "incident_energy_MeV" => system.incident_energy,
        "compound_A" => system.compound.A,
        "compound_Z" => system.compound.Z,
    )
end

"""
    symmetric_mass(system) -> Float64

The mass number of a symmetric split of the compound nucleus, `A₀ / 2`. Not generally an
integer; requested heavy-fragment ranges are validated against it.
"""
symmetric_mass(system::FissioningSystem) = system.compound.A / 2

const ELEMENT_SYMBOLS = (
    "n",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
)

"""
    element_symbol(Z) -> String

Chemical symbol for proton number `Z`, `0 ≤ Z ≤ 118`, with `"n"` at `Z = 0`. Any other `Z`
throws an `ArgumentError`: a system token built from a placeholder symbol would name no system.
"""
function element_symbol(Z::Integer)
    0 <= Z <= length(ELEMENT_SYMBOLS) - 1 ||
        throw(ArgumentError("no element symbol for Z = $Z; the table spans 0 to 118"))
    return ELEMENT_SYMBOLS[Z + 1]
end

function Base.show(io::IO, system::FissioningSystem)
    label = string(system.target.A, element_symbol(system.target.Z))
    if is_spontaneous(system)
        print(io, "FissioningSystem(", label, "(sf))")
    else
        print(
            io,
            "FissioningSystem(",
            label,
            "(",
            system.channel,
            ",f) at E = ",
            system.incident_energy,
            " MeV)",
        )
    end
    return nothing
end
