"""
    WahlReactionParameters

`Zₚ`-model parameters determined for one fission reaction by least squares, rather than estimated
from a trend across reactions. Fields are Table A of A. C. Wahl, *Atomic Data and Nuclear Data
Tables* **39**, 1 (1988), doi:10.1016/0092-640X(88)90016-2, p. 7, which fits the model to the
evaluated fractional independent yields of that reaction.

# Fields

- `assumed`: parameters the evaluation fixed by assumption, parenthesised in the table; they
  carry no uncertainty and are not measurements.
- `channels`: entrance channels the fit covers, so that a reaction reaching the same compound
  nucleus at a different excitation energy does not match.
- `shell`, `symmetric`: where the narrow dispersion `σ_Z(50)` applies and where the even-odd
  factors are unity, as footnotes a and c of Table A give them. Each covers the `Z = 50` closure
  and its complement by mass, so both legs appear where the table gives two.
"""
struct WahlReactionParameters
    label::String
    channels::Vector{String}
    Z_F::Int
    A_F::Int
    ΔZ140::Float64
    ΔZSL::Float64
    σZ::Float64
    σZ50::Float64
    F_Z::Float64
    F_N::Float64
    ΔA_Z::Float64
    ΔZmax::Float64
    shell::Vector{UnitRange{Int}}
    symmetric::UnitRange{Int}
    assumed::Vector{Symbol}

    # Explicit inner constructor taking exactly the concrete field types, with no conversion.
    function WahlReactionParameters(
        label::String,
        channels::Vector{String},
        Z_F::Int,
        A_F::Int,
        ΔZ140::Float64,
        ΔZSL::Float64,
        σZ::Float64,
        σZ50::Float64,
        F_Z::Float64,
        F_N::Float64,
        ΔA_Z::Float64,
        ΔZmax::Float64,
        shell::Vector{UnitRange{Int}},
        symmetric::UnitRange{Int},
        assumed::Vector{Symbol},
    )
        return new(
            label,
            channels,
            Z_F,
            A_F,
            ΔZ140,
            ΔZSL,
            σZ,
            σZ50,
            F_Z,
            F_N,
            ΔA_Z,
            ΔZmax,
            shell,
            symmetric,
            assumed,
        )
    end
end

function Base.show(io::IO, parameters::WahlReactionParameters)
    print(
        io,
        "WahlReactionParameters(",
        parameters.label,
        ", Z_F = ",
        parameters.Z_F,
        ", A_F = ",
        parameters.A_F,
        ")",
    )
    return nothing
end

"""
    WAHL_PER_REACTION

Table A of Wahl, At. Data Nucl. Data Tables 39, 1 (1988), doi:10.1016/0092-640X(88)90016-2,
keyed by the charge and mass of the fissioning nucleus, for the four reactions that evaluation
covers. Values with their uncertainties, as published:

| Parameter | U235T | U233T | PU239T | CF252S |
|---|---|---|---|---|
| `ΔZ(A' = 140)` | −0.511(5) | −0.519(8) | −0.544(9) | −0.420(20) |
| `∂ΔZ/∂A'` | −0.008(1) | −0.015(2) | −0.015(2) | −0.015(6) |
| `σ_Z` | 0.531(4) | 0.555(6) | 0.564(6) | 0.589(13) |
| `σ_Z(50)` | 0.33(2) | 0.36(3) | (0.35) | (0.35) |
| `F_Z` | 1.27(1) | 1.27(2) | 1.14(2) | 1.05(4) |
| `F_N` | 1.07(1) | 1.07(2) | 1.05(2) | (1.00) |
| `ΔA'_Z` | 0.9(2) | (0.9) | (0.9) | (0.7) |
| `ΔZ_max` | 0.7(1) | (0.7) | (0.7) | (0.5) |

Parenthesised values were assumed rather than determined and are listed in each entry's
`assumed` field. `σ_Z` is a single average over the mass range: the 1988 evaluation takes
`σ_Z(A) = σ_Z` outside the shell region and `σ_Z(50)` inside it. It is the model's Gaussian
width, which [`wahl_charge_distribution`](@ref) puts on the lattice footing of eq. (8) before
returning a [`ChargeDistribution`](@ref).
"""
const WAHL_PER_REACTION = Dict(
    (92, 236) => WahlReactionParameters(
        "U235T",
        ["nth", "nres"],
        92,
        236,
        -0.511,
        -0.008,
        0.531,
        0.33,
        1.27,
        1.07,
        0.9,
        0.7,
        [105:109, 125:129],
        105:129,
        Symbol[],
    ),
    (92, 234) => WahlReactionParameters(
        "U233T",
        ["nth", "nres"],
        92,
        234,
        -0.519,
        -0.015,
        0.555,
        0.36,
        1.27,
        1.07,
        0.9,
        0.7,
        [104:108, 123:128],
        104:128,
        [:ΔA_Z, :ΔZmax],
    ),
    (94, 240) => WahlReactionParameters(
        "PU239T",
        ["nth", "nres"],
        94,
        240,
        -0.544,
        -0.015,
        0.564,
        0.35,
        1.14,
        1.05,
        0.9,
        0.7,
        [109:113, 124:129],
        109:129,
        [:σZ50, :ΔA_Z, :ΔZmax],
    ),
    (98, 252) => WahlReactionParameters(
        "CF252S",
        ["sf"],
        98,
        252,
        -0.420,
        -0.015,
        0.589,
        0.35,
        1.05,
        1.00,
        0.7,
        0.5,
        [119:129],
        119:129,
        [:σZ50, :F_N, :ΔA_Z, :ΔZmax],
    ),
)

"""
    wahl_reaction_parameters(system) -> Union{WahlReactionParameters,Nothing}

The per-reaction `Zₚ` parameters for a fissioning system, or `nothing` where Wahl, At. Data
Nucl. Data Tables 39, 1 (1988), doi:10.1016/0092-640X(88)90016-2, covers no such reaction.

Matched on the fissioning nucleus and the entrance channel: ²⁴⁰Pu(sf) and ²³⁹Pu(n_th,f) reach
the same compound nucleus 6.53 MeV apart, over which eq. (17) of LA-13928 (2002),
doi:10.2172/809574, moves `σ_Z` by twelve times the uncertainty Table A quotes, so such a
reaction falls through to the systematics.
`nres` is admitted with `nth`, the two differing only by an incident energy that moves no
parameter at the sixth decimal.
"""
function wahl_reaction_parameters(system::FissioningSystem)
    parameters = get(WAHL_PER_REACTION, (system.compound.Z, system.compound.A), nothing)
    parameters === nothing && return nothing
    return system.channel in parameters.channels ? parameters : nothing
end

"""
    WahlSystematics(parameters, PE) -> WahlSystematics

Carry per-reaction parameters through the region functions of the `Zₚ` model. The regions and
their boundaries are those of LA-13928 (2002), doi:10.2172/809574; the parameters within them
are those of Table A of Wahl, At. Data Nucl. Data Tables 39, 1 (1988),
doi:10.1016/0092-640X(88)90016-2, for this reaction. The near-symmetry slope `SL50` and the wing
slopes, absent from Table A, come from the systematics.

`σ_Z` carries no slope, one average width applying outside the shell region, and the shell and
near-symmetry ranges of footnotes a and c take precedence over the boundaries of eq. (10). The
footnote ranges, given in fission-product mass, are applied to pre-neutron mass without a shift,
since `ν(A')` at the `Z = 50` closure, 0.2 to 0.7 neutrons, is below their one-mass-unit
resolution.
"""
function WahlSystematics(parameters::WahlReactionParameters, PE::Real)
    systematics = WahlSystematics(parameters.Z_F, parameters.A_F, PE)
    A_max = systematics.Bb
    B4 =
        (
            parameters.ΔZmax - parameters.ΔZ140 +
            A_max * systematics.SL50 +
            140 * parameters.ΔZSL
        ) / (systematics.SL50 + parameters.ΔZSL)
    return WahlSystematics(
        parameters.Z_F,
        parameters.A_F,
        Float64(PE),
        parameters.σZ,                 # σZ140: one average width, no slope
        parameters.ΔZ140,
        0.0,                           # σZSL
        parameters.ΔZSL,
        systematics.SL50,              # not in Table A
        parameters.σZ50,
        parameters.ΔZmax,
        systematics.σZSLW,             # not in Table A
        systematics.ΔZSLW,             # not in Table A
        parameters.F_Z,
        parameters.F_N,
        systematics.FZSLW,             # not in Table A
        systematics.FNSLW,             # not in Table A
        systematics.B2,
        B4,
        systematics.B5,
        systematics.B6,
        systematics.Ba,
        systematics.Bb,
        parameters.shell,
        parameters.symmetric,
    )
end

"""
    build_charge_distribution(table, system, masses; even_odd, charges) -> ChargeDistribution

The best isobaric charge distribution available for a fissioning system, resolved in three
layers, each falling through to the next where it cannot answer:

1. Per-reaction parameters. [`WAHL_PER_REACTION`](@ref) carries the `Zₚ` parameters of Wahl,
   At. Data Nucl. Data Tables 39, 1 (1988), doi:10.1016/0092-640X(88)90016-2, fitted by least
   squares to the measured fractional independent yields of four reactions.
2. The systematics. LA-13928 (2002), doi:10.2172/809574, eq. (17), for any nucleus inside
   [`WAHL_VALIDITY`](@ref); Fig. 18 of that report gives a reduced `χ²` of 7.9 against 2.9 for a
   per-reaction fit.
3. The conventional means `ΔZ = −0.5`, `σ_Z = 0.6` at every mass, the yield-weighted averages of
   the evaluated tables, without mass dependence.

The layer used is recorded with its citation in the returned distribution's `source`. An
evaluated table supplied through `charge_distribution_file` supersedes all three layers, and this
function is then not called.
"""
function build_charge_distribution(
    table::MassExcessTable,
    system::FissioningSystem,
    masses;
    even_odd::Bool = false,
    charges::Integer = 11,
)
    excitation = compound_nucleus_excitation(table, system)
    parameters = wahl_reaction_parameters(system)
    if parameters !== nothing && excitation !== nothing
        PE = value(excitation)
        if 0 <= PE <= WAHL_VALIDITY.excitation
            return wahl_charge_distribution(
                WahlSystematics(parameters, PE),
                collect(masses);
                even_odd = even_odd,
                charges = charges,
                provenance = "Zₚ model, per-reaction parameters for $(parameters.label), " *
                             "Wahl, At. Data Nucl. Data Tables 39, 1 (1988), Table A",
            )
        end
    end
    if is_wahl_applicable(table, system)
        distribution = wahl_charge_distribution(
            table,
            system,
            collect(masses);
            even_odd = even_odd,
            charges = charges,
        )
        distribution === nothing || return distribution
    end
    return mean_charge_distribution()
end
