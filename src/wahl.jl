"""
    WAHL_LOW_ENERGY

Coefficients of the `Zₚ` model of A. C. Wahl, *Systematics of Fission-Product Yields*, LA-13928
(2002), doi:10.2172/809574, Table 2, for the low-energy branch `PE ≤ 8 MeV`. Each entry is
`(P₁ … P₅)` of eq. (17),

```
Par = P₁ + P₂[Z_F − 92] + P₃[A_F − 236] + P₄[PE − 6.551] + P₅[A_F − 236]²
```

The peak-region slope `FZ SL = 0.0030` of Table 2 is not
carried, since eq. (11d) makes `F_Z(A')` constant there; the `F_N` slope of footnote b,
`−0.0006(10)`, is consistent with zero.
"""
const WAHL_LOW_ENERGY = (
    σZ140 = (0.566, 0.0, 0.0064, 0.0109, 0.0),
    ΔZ140 = (-0.487, 0.0, 0.0180, 0.0, -0.00203),
    σZSL = (-0.0038, 0.0, 0.0, 0.0, 0.0),
    ΔZSL = (-0.0080, 0.0, 0.0, 0.0, 0.0),
    SL50 = (0.191, 0.0, -0.0076, 0.0, 0.0),
    σZ50 = (0.356, 0.060, 0.0, 0.0, 0.0),
    ΔZmax = (0.699, 0.0, 0.0, 0.0, 0.0),
    σZSLW = (-0.045, 0.0094, 0.0, 0.0, 0.0),
    ΔZSLW = (0.0, -0.0045, 0.0, 0.0, 0.0),
    FZ140 = (1.207, 0.0, -0.0420, 0.0, 0.0022),
    FN140 = (1.076, 0.0, 0.0, 0.0, 0.0),
    FZSLW = (0.159, -0.028, 0.0, 0.0, 0.0),
    FNSLW = (0.039, 0.0, 0.0, 0.0, 0.0),
)

"""
    WAHL_VALIDITY

The range of fissioning nuclides the systematics was fitted over: `Z_F ∈ [90, 98]`,
`A_F ∈ [230, 252]`, and, for the low-energy branch used here, `PE ≤ 8 MeV`.
"""
const WAHL_VALIDITY = (charge = 90:98, mass = 230:252, excitation = 8.0)

"""
    WahlSystematics <: ZpModel

The `Zₚ` model evaluated for one fissioning nucleus from the systematics across reactions: the
parameters of eq. (17) of LA-13928 (2002), doi:10.2172/809574, and the mass-number boundaries of
eq. (10) that divide `A'` into wing, peak and near-symmetry regions. `A'` is the precursor mass,
the mass before prompt-neutron emission, and so the index of primary fragments. `ΔZ` and `σ_Z`
are quoted for the heavy side, the convention of eqs. (9c) and (9d).

`Ba` and `Bb` bound the interval where `Zₚ` crosses the `Z = 50` shell closure, with
`Ba = A_F − A'_max` and `Bb = A'_max` as in Figs. 18b, 19b and 20b; the printed eq. (10) swaps
the two, and only the assignment of the figures is continuous with the adjacent peak regions.
Where a reaction has its own parameters, [`Wahl1988`](@ref) takes precedence;
[`charge_model`](@ref) selects between the two.
"""
struct WahlSystematics <: ZpModel
    Z_F::Int
    A_F::Int
    PE::Float64
    σZ140::Float64
    ΔZ140::Float64
    σZSL::Float64
    ΔZSL::Float64
    SL50::Float64
    σZ50::Float64
    ΔZmax::Float64
    σZSLW::Float64
    ΔZSLW::Float64
    FZ140::Float64
    FN140::Float64
    FZSLW::Float64
    FNSLW::Float64
    B2::Float64
    B4::Float64
    B5::Float64
    B6::Float64
    Ba::Float64
    Bb::Float64
    neutron_pairing::Bool

    # Explicit inner constructor converting every argument to its concrete field type.
    function WahlSystematics(
        Z_F::Integer,
        A_F::Integer,
        PE::Real,
        σZ140::Real,
        ΔZ140::Real,
        σZSL::Real,
        ΔZSL::Real,
        SL50::Real,
        σZ50::Real,
        ΔZmax::Real,
        σZSLW::Real,
        ΔZSLW::Real,
        FZ140::Real,
        FN140::Real,
        FZSLW::Real,
        FNSLW::Real,
        B2::Real,
        B4::Real,
        B5::Real,
        B6::Real,
        Ba::Real,
        Bb::Real,
        neutron_pairing::Bool,
    )
        return new(
            Int(Z_F),
            Int(A_F),
            Float64(PE),
            Float64(σZ140),
            Float64(ΔZ140),
            Float64(σZSL),
            Float64(ΔZSL),
            Float64(SL50),
            Float64(σZ50),
            Float64(ΔZmax),
            Float64(σZSLW),
            Float64(ΔZSLW),
            Float64(FZ140),
            Float64(FN140),
            Float64(FZSLW),
            Float64(FNSLW),
            Float64(B2),
            Float64(B4),
            Float64(B5),
            Float64(B6),
            Float64(Ba),
            Float64(Bb),
            neutron_pairing,
        )
    end
end

_wahl_parameter(coefficients, Z_F::Integer, A_F::Integer, PE::Real) =
    coefficients[1] +
    coefficients[2] * (Z_F - 92) +
    coefficients[3] * (A_F - 236) +
    coefficients[4] * (PE - 6.551) +
    coefficients[5] * (A_F - 236)^2

"""
    WahlSystematics(Z_F, A_F, PE; neutron_pairing = false) -> WahlSystematics

Evaluate the `Zₚ` model for a fissioning nucleus `(Z_F, A_F)` at precursor excitation energy
`PE` in MeV — zero for spontaneous fission, `Sₙ + Eₙ` for neutron-induced fission. Throws for
a nucleus or an excitation energy outside [`WAHL_VALIDITY`](@ref), the range the systematics was
fitted over. `neutron_pairing` is the fragment-level convention of
[`fragment_charge_yields`](@ref).
"""
function WahlSystematics(
    Z_F::Integer,
    A_F::Integer,
    PE::Real;
    neutron_pairing::Bool = false,
)
    Z_F in WAHL_VALIDITY.charge || throw(
        ArgumentError("the Zₚ systematics is fitted for Z_F in $(WAHL_VALIDITY.charge), \
                       got $Z_F"),
    )
    A_F in WAHL_VALIDITY.mass || throw(
        ArgumentError("the Zₚ systematics is fitted for A_F in $(WAHL_VALIDITY.mass), \
                       got $A_F"),
    )
    (isfinite(PE) && PE >= 0) || throw(
        ArgumentError("the precursor excitation energy must be finite and non-negative, \
                             got $PE MeV"),
    )
    PE <= WAHL_VALIDITY.excitation ||
        throw(ArgumentError("the low-energy branch of the Zₚ systematics holds to \
                       PE = $(WAHL_VALIDITY.excitation) MeV, got $PE MeV"))

    evaluate(key) = _wahl_parameter(getproperty(WAHL_LOW_ENERGY, key), Z_F, A_F, PE)
    σZ140 = evaluate(:σZ140)
    ΔZ140 = evaluate(:ΔZ140)
    σZSL = evaluate(:σZSL)
    ΔZSL = evaluate(:ΔZSL)
    SL50 = evaluate(:SL50)
    σZ50 = evaluate(:σZ50)
    ΔZmax = evaluate(:ΔZmax)
    σZSLW = evaluate(:σZSLW)
    ΔZSLW = evaluate(:ΔZSLW)
    FZ140 = evaluate(:FZ140)
    FN140 = evaluate(:FN140)
    FZSLW = evaluate(:FZSLW)
    FNSLW = evaluate(:FNSLW)

    # Eq. (12): the A' at which Zₚ reaches 50, a blend of two estimates AK1 and AK2 weighted by
    # F₁. The blend places `Bb`, the mass of maximum ΔZ, closer to the evaluated tables than AK2.
    F₁ = clamp((250 - A_F) / 14, 0.0, 1.0)
    A_max = F₁ * (50 * (A_F / Z_F) - ΔZmax / SL50) + (1 - F₁) * ((50 - ΔZmax) * (A_F / Z_F))

    B2 = 77 + 0.036 * (A_F - 236)
    B4 = (ΔZmax - ΔZ140 + A_max * SL50 + 140 * ΔZSL) / (SL50 + ΔZSL)
    return WahlSystematics(
        Int(Z_F),
        Int(A_F),
        Float64(PE),
        σZ140,
        ΔZ140,
        σZSL,
        ΔZSL,
        SL50,
        σZ50,
        ΔZmax,
        σZSLW,
        ΔZSLW,
        FZ140,
        FN140,
        FZSLW,
        FNSLW,
        B2,
        B4,
        A_F - B2,
        A_F - 70.0,
        A_F - A_max,
        A_max,
        neutron_pairing,
    )
end

"""
    is_wahl_applicable(table, system) -> Bool

Whether the `Zₚ` systematics can be evaluated for a fissioning system at all: the compound
nucleus inside the fitted `Z_F` and `A_F` ranges of [`WAHL_VALIDITY`](@ref), and a precursor
excitation energy the mass table can form and the low-energy branch covers. Queried before
construction, so that a caller can omit the systematics where it does not apply instead of
catching the constructor's exception.
"""
function is_wahl_applicable(table::MassExcessTable, system::FissioningSystem)
    system.compound.Z in WAHL_VALIDITY.charge || return false
    system.compound.A in WAHL_VALIDITY.mass || return false
    excitation = compound_nucleus_excitation(table, system)
    excitation === nothing && return false
    PE = value(excitation)
    return isfinite(PE) && 0 <= PE <= WAHL_VALIDITY.excitation
end

"""
    WahlSystematics(table, system; neutron_pairing = false) -> Union{WahlSystematics,Nothing}

The systematics for a fissioning system, taking the precursor excitation energy from the
compound nucleus: zero for spontaneous fission, `Sₙ + Eₙ` otherwise. Returns `nothing` when the
mass table cannot form `Sₙ`.
"""
function WahlSystematics(
    table::MassExcessTable,
    system::FissioningSystem;
    neutron_pairing::Bool = false,
)
    excitation = compound_nucleus_excitation(table, system)
    excitation === nothing && return nothing
    return WahlSystematics(
        system.compound.Z,
        system.compound.A,
        value(excitation);
        neutron_pairing = neutron_pairing,
    )
end

_wahl_peak_polarization(w::WahlSystematics, A) = w.ΔZ140 + w.ΔZSL * (A - 140)
_wahl_peak_dispersion(w::WahlSystematics, A) = w.σZ140 + w.σZSL * (A - 140)

"""
    wahl_polarization(systematics, A) -> Float64

Charge distribution `ΔZ` of the heavy fragment at pre-neutron mass `A`, eqs. (11a), (12b),
(12d), (12f), (13b) and (14b). Defined for `A ≥ A_F/2`; the light partner takes `−ΔZ` at its
complement, as eq. (9d) states and [`most_probable_charge`](@ref) applies.
"""
function wahl_polarization(systematics::WahlSystematics, A::Real)
    _wahl_check_heavy(systematics, A)
    A >= systematics.B6 && return _wahl_peak_polarization(systematics, systematics.B5)
    A >= systematics.B5 && return _wahl_peak_polarization(systematics, systematics.B5) -
           systematics.ΔZSLW * (A - systematics.B5)
    A >= systematics.B4 && return _wahl_peak_polarization(systematics, A)
    edge =
        _wahl_peak_polarization(systematics, systematics.B4) +
        systematics.SL50 * (systematics.B4 - systematics.Bb)
    A >= systematics.Bb && return _wahl_peak_polarization(systematics, systematics.B4) +
           systematics.SL50 * (systematics.B4 - A)
    # Ba to Bb: the Z = 50 crossing, a straight line between values equal and opposite by
    # complementarity, so ΔZ vanishes at symmetry as charge conservation requires.
    return -edge + 2 * edge * (A - systematics.Ba) / (systematics.Bb - systematics.Ba)
end

"""
    wahl_dispersion(systematics, A) -> Float64

Dispersion `σ_Z` of the isobaric charge distribution at pre-neutron mass `A`, eqs. (11b), (12c),
(12e), (12g), (13d) and (14c). Defined for `A ≥ A_F/2`.

This is the model parameter; the second moment of the yields over the integer charge lattice is
`√(σ_Z² + 1/12)`, eq. (8) of Wahl, At. Data Nucl. Data Tables 39, 1 (1988),
doi:10.1016/0092-640X(88)90016-2. Unlike `ΔZ`, `σ_Z` is a step function: eqs. (12c) and (12g) hold
`σ_Z(50)` flat over the two legs adjoining the `Z = 50` crossing while the crossing keeps the
peak value, and eq. (14c) returns the far wing to `σ_Z(B5)` (Fig. 18a, p. 31).
"""
function wahl_dispersion(systematics::WahlSystematics, A::Real)
    _wahl_check_heavy(systematics, A)
    A >= systematics.B6 && return _wahl_peak_dispersion(systematics, systematics.B5)
    A >= systematics.B5 && return _wahl_peak_dispersion(systematics, systematics.B5) +
           systematics.σZSLW * (A - systematics.B5)
    A >= systematics.B4 && return _wahl_peak_dispersion(systematics, A)
    A >= systematics.Bb && return systematics.σZ50
    return _wahl_peak_dispersion(systematics, systematics.Bb)
end

"""
    wahl_even_odd_factors(systematics, A) -> (F_Z, F_N)

The even-odd proton and neutron factors at heavy-fragment mass `A`, eqs. (11d), (11e), (12a),
(13f), (13h), (14d) and (14e).

Both are `1` through the near-symmetry region, where eq. (12a) sets `F(A) = 1`, take their peak
values through the peak and far-wing regions, and carry a slope only in the wing regions.
"""
function wahl_even_odd_factors(systematics::WahlSystematics, A::Real)
    _wahl_check_heavy(systematics, A)
    A >= systematics.B6 && return (systematics.FZ140, systematics.FN140)
    A >= systematics.B5 && return (
        systematics.FZ140 + systematics.FZSLW * (A - systematics.B5),
        systematics.FN140 + systematics.FNSLW * (A - systematics.B5),
    )
    A >= systematics.B4 && return (systematics.FZ140, systematics.FN140)
    return (1.0, 1.0)
end

function _wahl_check_heavy(systematics::WahlSystematics, A::Real)
    2 * A >= systematics.A_F ||
        throw(ArgumentError("the Zₚ systematics is evaluated for the heavy fragment, \
                       A ≥ A_F/2 = $(systematics.A_F / 2), got $A"))
    return nothing
end

_heavy_parameters(systematics::WahlSystematics, A′::Real) = (
    wahl_polarization(systematics, A′),
    wahl_dispersion(systematics, A′),
    wahl_even_odd_factors(systematics, A′)...,
)

function Base.show(io::IO, systematics::WahlSystematics)
    print(
        io,
        "WahlSystematics(Z_F = ",
        systematics.Z_F,
        ", A_F = ",
        systematics.A_F,
        ", PE = ",
        round(systematics.PE; digits = 3),
        " MeV)",
    )
    return nothing
end
