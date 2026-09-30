"""
    Wahl1988 <: ZpModel

The `Zₚ` model of A. C. Wahl, *At. Data Nucl. Data Tables* **39**, 1 (1988),
doi:10.1016/0092-640X(88)90016-2, with the parameters that evaluation determined for one fission
reaction by least squares (Table A). The functions of `A'`:

- `ΔZ(A'_H)` is the straight line `ΔZ(140) + ∂ΔZ/∂A' (A'_H − 140)` down to the junction `A'_J`.
  Below it a steep branch rises across the `Zₚ = 50` line, and below `A'_m`, where it reaches
  `ΔZ_max`, a straight line falls to zero at `A_F/2` (Fig. 2). For CF252S the steep branch
  reaches `ΔZ_max` only below `A_F/2`, so it runs down to symmetry and `ΔZ(A_F/2) ≈ 0.49`; the
  heavy and the light product of mass `A_F/2` then centre on different charges, as Table IV has
  them. The steep branch passes through
  the point where the `Zₚ = 50` line meets `ΔZ = 0`, `A'_{50} = 50 A_F/Z_F`. It also passes
  through the point `ΔA'_Z` above `A'_P` at the height of `P`, where `P` is the intersection of the
  `Zₚ = 50` line with the extrapolated peak line; that point is the junction.
- `σ_Z = σ̄₅₀` on the steep branch, `A'_m ≤ A'_H ≤ A'_J`, and `σ̄_Z` elsewhere.
- `F(A) = 1` for `A'_H ≤ A'_J`, and `F̄_Z`, `F̄_N` above.

The figure and the text leave the direction of the `ΔA'_Z` displacement and the region
boundaries to be read. The reading here is the one that reproduces the evaluation's own
calculations: every `Zₚ(A)` of Tables I–IV to 0.0007 from the tabulated `ν̄_A`, and every
calculated `Z̄(A)` and `RMS(A)` to 0.0018. Footnote a of Table A lists the `σ̄₅₀` masses as
products. The steep-branch rule reproduces those lists, except that for U235T it excludes
`A = 109`, as the tables themselves do. See [`WAHL_1988`](@ref) for the parameter values.
"""
struct Wahl1988 <: ZpModel
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
    A_50::Float64
    A_J::Float64
    A_m::Float64
    steep::Float64
    neutron_pairing::Bool

    function Wahl1988(
        label::AbstractString,
        channels::Vector{String},
        Z_F::Integer,
        A_F::Integer,
        ΔZ140::Real,
        ΔZSL::Real,
        σZ::Real,
        σZ50::Real,
        F_Z::Real,
        F_N::Real,
        ΔA_Z::Real,
        ΔZmax::Real;
        neutron_pairing::Bool = false,
    )
        (σZ > 0 && σZ50 > 0) || throw(ArgumentError("the widths must be positive"))
        (F_Z > 0 && F_N > 0) ||
            throw(ArgumentError("the even-odd factors must be positive"))
        k = Z_F / A_F
        # P: the Zₚ = 50 line, ΔZ = 50 − kA', meets the extrapolated peak line.
        A_P = (50 - ΔZ140 + 140 * ΔZSL) / (k + ΔZSL)
        ΔZ_P = 50 - k * A_P
        A_50 = 50 / k
        A_J = A_P + ΔA_Z
        steep = -ΔZ_P / (A_50 - A_J)
        A_m = A_50 + ΔZmax / steep
        A_m < A_J || throw(
            ArgumentError(
                "the steep branch must lie below the junction: A'_m = $A_m, A'_J = $A_J",
            ),
        )
        return new(
            String(label),
            channels,
            Int(Z_F),
            Int(A_F),
            Float64(ΔZ140),
            Float64(ΔZSL),
            Float64(σZ),
            Float64(σZ50),
            Float64(F_Z),
            Float64(F_N),
            Float64(ΔA_Z),
            Float64(ΔZmax),
            A_50,
            A_J,
            A_m,
            steep,
            neutron_pairing,
        )
    end
end

function _heavy_parameters(model::Wahl1988, A′::Real)
    ΔZ = if A′ >= model.A_J
        model.ΔZ140 + model.ΔZSL * (A′ - 140)
    elseif A′ >= model.A_m
        model.steep * (A′ - model.A_50)
    else
        model.ΔZmax * (A′ - model.A_F / 2) / (model.A_m - model.A_F / 2)
    end
    σ_Z = model.A_m <= A′ <= model.A_J ? model.σZ50 : model.σZ
    F_Z, F_N = A′ <= model.A_J ? (1.0, 1.0) : (model.F_Z, model.F_N)
    return ΔZ, σ_Z, F_Z, F_N
end

function Base.show(io::IO, model::Wahl1988)
    print(io, "Wahl1988(", model.label, ", Z_F = ", model.Z_F, ", A_F = ", model.A_F, ")")
    return nothing
end

"""
    WAHL_1988_TABLE_A

Table A of Wahl, *At. Data Nucl. Data Tables* **39**, 1 (1988), doi:10.1016/0092-640X(88)90016-2,
p. 7, as printed: for each reaction, each parameter as `(value, uncertainty)`, with `nothing` in
place of the uncertainty where the value was assumed rather than fitted (parenthesised in the
table).

| Parameter | U235T | U233T | PU239T | CF252S |
|---|---|---|---|---|
| `ΔZ(A' = 140)` | −0.511(5) | −0.519(8) | −0.544(9) | −0.420(20) |
| `∂ΔZ/∂A'` | −0.008(1) | −0.015(2) | −0.015(2) | −0.015(6) |
| `σ̄_Z` | 0.531(4) | 0.555(6) | 0.564(6) | 0.589(13) |
| `σ̄₅₀` | 0.33(2) | 0.36(3) | (0.35) | (0.35) |
| `F̄_Z` | 1.27(1) | 1.27(2) | 1.14(2) | 1.05(4) |
| `F̄_N` | 1.07(1) | 1.07(2) | 1.05(2) | (1.00) |
| `ΔA'_Z` | 0.9(2) | (0.9) | (0.9) | (0.7) |
| `ΔZ_max` | 0.7(1) | (0.7) | (0.7) | (0.5) |

The values are rounded to the uncertainty. [`WAHL_1988`](@ref) carries the values the tables were
calculated with.
"""
const WAHL_1988_TABLE_A = Dict(
    "U235T" => (
        ΔZ140 = (-0.511, 0.005),
        ΔZSL = (-0.008, 0.001),
        σZ = (0.531, 0.004),
        σZ50 = (0.33, 0.02),
        F_Z = (1.27, 0.01),
        F_N = (1.07, 0.01),
        ΔA_Z = (0.9, 0.2),
        ΔZmax = (0.7, 0.1),
    ),
    "U233T" => (
        ΔZ140 = (-0.519, 0.008),
        ΔZSL = (-0.015, 0.002),
        σZ = (0.555, 0.006),
        σZ50 = (0.36, 0.03),
        F_Z = (1.27, 0.02),
        F_N = (1.07, 0.02),
        ΔA_Z = (0.9, nothing),
        ΔZmax = (0.7, nothing),
    ),
    "PU239T" => (
        ΔZ140 = (-0.544, 0.009),
        ΔZSL = (-0.015, 0.002),
        σZ = (0.564, 0.006),
        σZ50 = (0.35, nothing),
        F_Z = (1.14, 0.02),
        F_N = (1.05, 0.02),
        ΔA_Z = (0.9, nothing),
        ΔZmax = (0.7, nothing),
    ),
    "CF252S" => (
        ΔZ140 = (-0.420, 0.020),
        ΔZSL = (-0.015, 0.006),
        σZ = (0.589, 0.013),
        σZ50 = (0.35, nothing),
        F_Z = (1.05, 0.04),
        F_N = (1.00, nothing),
        ΔA_Z = (0.7, nothing),
        ΔZmax = (0.5, nothing),
    ),
)

"""
    WAHL_1988

The four reactions of Wahl, *At. Data Nucl. Data Tables* **39**, 1 (1988),
doi:10.1016/0092-640X(88)90016-2, as [`Wahl1988`](@ref) models keyed by the charge and mass of
the fissioning nucleus.

The parameters are those the evaluation's Tables I–IV were calculated with. They are recovered by
least squares from the tabulated `Zₚ(A)`, `Z̄(A)` and `RMS(A)`, and each lies within the rounding
of the value [`WAHL_1988_TABLE_A`](@ref) prints:

- `∂ΔZ/∂A'` is −0.0079 for U235T, and −0.0146 for U233T and PU239T;
- `ΔA'_Z` and `ΔZ_max` are 0.88 and 0.68 for U235T;
- `F̄_Z` and `F̄_N` are 1.144 and 1.046 for PU239T.

The rest are as printed. With these values the tables are reproduced at their printed precision;
with the printed ones, `Zₚ` departs by up to 0.013 at the extreme masses.

`channels` names the entrance channels each fit covers; `nres` shares the thermal fit.
"""
const WAHL_1988 = Dict(
    (92, 236) => Wahl1988(
        "U235T",
        ["nth", "nres"],
        92,
        236,
        -0.511,
        -0.0079,
        0.531,
        0.33,
        1.27,
        1.07,
        0.88,
        0.68,
    ),
    (92, 234) => Wahl1988(
        "U233T",
        ["nth", "nres"],
        92,
        234,
        -0.519,
        -0.0146,
        0.555,
        0.36,
        1.27,
        1.07,
        0.9,
        0.7,
    ),
    (94, 240) => Wahl1988(
        "PU239T",
        ["nth", "nres"],
        94,
        240,
        -0.544,
        -0.0146,
        0.564,
        0.35,
        1.144,
        1.046,
        0.9,
        0.7,
    ),
    (98, 252) => Wahl1988(
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
    ),
)

"""
    Wahl1988(label::AbstractString; printed = false, neutron_pairing = false) -> Wahl1988

The model of one reaction by its label, `"U235T"`, `"U233T"`, `"PU239T"` or `"CF252S"`. The
values are those of [`WAHL_1988`](@ref), or with `printed = true` those of
[`WAHL_1988_TABLE_A`](@ref).
"""
function Wahl1988(
    label::AbstractString;
    printed::Bool = false,
    neutron_pairing::Bool = false,
)
    index = findfirst(m -> m.label == label, collect(values(WAHL_1988)))
    index === nothing && throw(
        ArgumentError(
            "no reaction $(repr(label)) in Wahl (1988); it covers U235T, U233T, PU239T and CF252S",
        ),
    )
    m = collect(values(WAHL_1988))[index]
    p =
        printed ? map(first, WAHL_1988_TABLE_A[m.label]) :
        (
            ΔZ140 = m.ΔZ140,
            ΔZSL = m.ΔZSL,
            σZ = m.σZ,
            σZ50 = m.σZ50,
            F_Z = m.F_Z,
            F_N = m.F_N,
            ΔA_Z = m.ΔA_Z,
            ΔZmax = m.ΔZmax,
        )
    return Wahl1988(
        m.label,
        m.channels,
        m.Z_F,
        m.A_F,
        p.ΔZ140,
        p.ΔZSL,
        p.σZ,
        p.σZ50,
        p.F_Z,
        p.F_N,
        p.ΔA_Z,
        p.ΔZmax;
        neutron_pairing = neutron_pairing,
    )
end

"""
    wahl_1988(system; neutron_pairing = false) -> Union{Wahl1988,Nothing}

The per-reaction model for a fissioning system, or `nothing` where Wahl (1988) covers no such
reaction. Matched on the fissioning nucleus and the entrance channel: ²⁴⁰Pu(sf) and
²³⁹Pu(n_th,f) reach the same compound nucleus 6.53 MeV apart, and eq. (17) of LA-13928 (2002),
doi:10.2172/809574, moves `σ_Z` between them by twelve times the uncertainty Table A quotes, so
such a reaction takes the systematics instead.
"""
function wahl_1988(system::FissioningSystem; neutron_pairing::Bool = false)
    model = get(WAHL_1988, (system.compound.Z, system.compound.A), nothing)
    (model === nothing || !(system.channel in model.channels)) && return nothing
    return Wahl1988(model.label; neutron_pairing = neutron_pairing)
end

"""
    charge_model(table, system; neutron_pairing = false) -> ChargeModel

The best charge model available for a fissioning system without an evaluated table, resolved in
three layers:

1. [`Wahl1988`](@ref), the reaction's own least-squares parameters, for the four reactions of
   Wahl (1988);
2. [`WahlSystematics`](@ref), eq. (17) of LA-13928 (2002), for any nucleus inside
   [`WAHL_VALIDITY`](@ref); Fig. 18 of that report gives a reduced `χ²` of 7.9 for it against
   2.9 for a per-reaction fit;
3. [`mean_charge_distribution`](@ref), `ΔZ = −0.5` and `σ_Z = 0.6` at every mass.

`neutron_pairing` is passed to the `Zₚ` models; see [`fragment_charge_yields`](@ref).
"""
function charge_model(
    table::MassExcessTable,
    system::FissioningSystem;
    neutron_pairing::Bool = false,
)
    reaction = wahl_1988(system; neutron_pairing = neutron_pairing)
    reaction === nothing || return reaction
    is_wahl_applicable(table, system) || return mean_charge_distribution()
    systematics = WahlSystematics(table, system; neutron_pairing = neutron_pairing)
    return systematics === nothing ? mean_charge_distribution() : systematics
end
