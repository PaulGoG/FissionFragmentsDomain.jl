"""
    ZpModel

The `Zₚ` model of A. C. Wahl: the dispersion of fission-product yields in `Z` at each mass,
a Gaussian integrated over unit charge intervals and modulated by even-odd proton and neutron
factors. Eq. (7) of Wahl, *At. Data Nucl. Data Tables* **39**, 1 (1988),
doi:10.1016/0092-640X(88)90016-2, and eq. (9) of LA-13928 (2002), doi:10.2172/809574:

```
FI(A,Z) = ½ F(A) N(A) [erf(V) − erf(W)],   V, W = (Z − Zₚ ± ½) / (σ_Z √2)
```

Every parameter is a function of the precursor mass `A'`, the mass before prompt-neutron
emission, and is defined for the heavy side, `A' ≥ A_F/2`. The light side follows by
complementarity, eq. (7d): `Zₚ(A'_L) = A'_L Z_F/A_F − ΔZ(A_F − A'_L)`, with `σ_Z` and the
even-odd factors taken at the complement. `N(A)` renormalises the yields at each mass to unit sum,
which the even-odd factors would otherwise break.

The two members are [`Wahl1988`](@ref), the per-reaction evaluation, and
[`WahlSystematics`](@ref), the systematics across reactions.

Wahl's yields are yields of products, with the neutron factor on the parity of the product's
`N`. A primary fragment's `N` differs from its product's by the evaporated multiplicity, so at
fragment level [`fragment_charge_yields`](@ref) applies the proton factor alone unless the model
was built with `neutron_pairing = true`. [`fractional_independent_yields`](@ref) evaluates the
model as Wahl defined it, for products.
"""
abstract type ZpModel <: ChargeModel end

# Subtypes provide `_heavy_parameters(model, A′H) -> (ΔZ, σ_Z, F_Z, F_N)` for A′H ≥ A_F/2, and
# the fields Z_F, A_F and neutron_pairing.

_heavy_mass(model::ZpModel, A′::Real) = 2 * A′ >= model.A_F ? A′ : model.A_F - A′

"""
    most_probable_charge(model::ZpModel, A′) -> Float64

`Zₚ` at precursor mass `A′`, heavy or light: eqs. (7c) and (7d) of Wahl (1988).
"""
function most_probable_charge(model::ZpModel, A′::Real)
    ΔZ = first(_heavy_parameters(model, _heavy_mass(model, A′)))
    UCD = A′ * model.Z_F / model.A_F
    return 2 * A′ >= model.A_F ? UCD + ΔZ : UCD - ΔZ
end

"""
    charge_polarization(model::ZpModel, A′) -> Float64

`ΔZ(A′)` for a heavy precursor mass `A′ ≥ A_F/2`, the convention of
[`ChargeDistribution`](@ref).
"""
function charge_polarization(model::ZpModel, A′::Real)
    2 * A′ >= model.A_F || throw(
        ArgumentError(
            "ΔZ is quoted for the heavy side, A′ ≥ A_F/2 = $(model.A_F / 2); got $A′",
        ),
    )
    return first(_heavy_parameters(model, A′))
end

"""
    charge_dispersion(model::ZpModel, A′) -> Float64

The width parameter `σ_Z` at precursor mass `A′`, heavy or light. It is the Gaussian's own
width; the second moment of the yields over the charge lattice is `√(σ_Z² + 1/12)`, eq. (8) of
Wahl (1988).
"""
charge_dispersion(model::ZpModel, A′::Real) =
    _heavy_parameters(model, _heavy_mass(model, A′))[2]

"""
    even_odd_factors(model::ZpModel, A′) -> (F_Z, F_N)

The even-odd proton and neutron factors at precursor mass `A′`, heavy or light; `(1, 1)` where
the model sets `F(A) = 1` near symmetry.
"""
function even_odd_factors(model::ZpModel, A′::Real)
    _, _, F_Z, F_N = _heavy_parameters(model, _heavy_mass(model, A′))
    return (F_Z, F_N)
end

function _pairing_factor(Z::Integer, N::Integer, F_Z::Real, F_N::Real)
    iseven(Z) && iseven(N) && return F_Z * F_N
    iseven(Z) && return F_Z / F_N
    iseven(N) && return F_N / F_Z
    return 1 / (F_Z * F_N)
end

# Eq. (7) over every charge within eight units of Zₚ, where σ_Z < 0.8 leaves no weight beyond,
# renormalised to unit sum. `A` sets the neutron parity; `F_N = 1` removes the neutron factor.
function _lattice_yields(Zₚ::Real, σ_Z::Real, F_Z::Real, F_N::Real, A::Integer)
    σ_Z > 0 || throw(ArgumentError("the charge dispersion must be positive, got $σ_Z"))
    charges = (floor(Int, Zₚ) - 8):(ceil(Int, Zₚ) + 8)
    scale = σ_Z * sqrt(2)
    yields = [
        0.5 *
        _pairing_factor(Z, A - Z, F_Z, F_N) *
        (erf((Z - Zₚ + 0.5) / scale) - erf((Z - Zₚ - 0.5) / scale)) for Z in charges
    ]
    yields ./= sum(yields)
    return charges, yields
end

"""
    fractional_independent_yields(model, A, ν̄) -> (charges, FI)

The fractional independent yields of the products of mass `A`, formed from precursors of mass
`A′ = A + ν̄` (eq. (2) of Wahl, 1988), exactly as Wahl defines them: both even-odd factors, on the
parities of the product's `Z` and `N = A − Z`. This is the form the calculated columns of
Tables I–IV of that evaluation tabulate. `ν̄ = 0` evaluates the model at `A` itself.

# Examples

```jldoctest
julia> charges, FI = fractional_independent_yields(WAHL_1988[(98, 252)], 140, 1.733);

julia> round.(charge_moments(charges, FI); digits = 3)   # Table IV prints 54.660 and 0.665
(54.659, 0.665)
```
"""
function fractional_independent_yields(model::ZpModel, A::Integer, ν̄::Real)
    (isfinite(ν̄) && ν̄ >= 0) ||
        throw(ArgumentError("ν̄ must be finite and non-negative, got $ν̄"))
    A′ = A + ν̄
    ΔZ, σ_Z, F_Z, F_N = _heavy_parameters(model, _heavy_mass(model, A′))
    UCD = A′ * model.Z_F / model.A_F
    Zₚ = 2 * A′ >= model.A_F ? UCD + ΔZ : UCD - ΔZ
    return _lattice_yields(Zₚ, σ_Z, F_Z, F_N, A)
end

"""
    fragment_charge_yields(model, A) -> (charges, p)

The charge distribution of primary fragments of mass `A`, before prompt-neutron emission: the
model evaluated at `A′ = A`. The proton factor applies on the fragment's `Z`, which neutron
emission leaves unchanged. The neutron factor applies on the fragment's `N` only for a model
built with `neutron_pairing = true`, since Wahl fitted it to the parity of the product's `N`.
The yields sum to one over the returned charges.
"""
function fragment_charge_yields(model::ZpModel, A::Integer)
    ΔZ, σ_Z, F_Z, F_N = _heavy_parameters(model, _heavy_mass(model, A))
    UCD = A * model.Z_F / model.A_F
    Zₚ = 2 * A >= model.A_F ? UCD + ΔZ : UCD - ΔZ
    return _lattice_yields(Zₚ, σ_Z, F_Z, model.neutron_pairing ? F_N : 1.0, A)
end

"""
    charge_moments(charges, yields) -> (Z̄, RMS)

The mean charge and the root-mean-square dispersion of a charge distribution normalised to unit
sum, eqs. (8a) and (8b) of Wahl (1988).
"""
function charge_moments(charges::AbstractVector{<:Integer}, yields::AbstractVector{<:Real})
    length(charges) == length(yields) ||
        throw(DimensionMismatch("charges and yields differ in length"))
    Z̄ = sum(Z * y for (Z, y) in zip(charges, yields))
    return Z̄, sqrt(sum((Z - Z̄)^2 * y for (Z, y) in zip(charges, yields)))
end

function _check_system(model::ZpModel, system::FissioningSystem)
    (system.compound.Z == model.Z_F && system.compound.A == model.A_F) || throw(
        ArgumentError(
            "the Zₚ model is for Z_F = $(model.Z_F), A_F = $(model.A_F), not for " *
            "$(system.compound)",
        ),
    )
    return nothing
end

fragment_most_probable_charge(model::ZpModel, system::FissioningSystem, A::Integer) =
    (_check_system(model, system); most_probable_charge(model, A))

function fragment_charge_probability(
    model::ZpModel,
    system::FissioningSystem,
    fragment::Nuclide,
)
    _check_system(model, system)
    charges, yields = fragment_charge_yields(model, fragment.A)
    index = fragment.Z - first(charges) + 1
    return 1 <= index <= length(yields) ? yields[index] : 0.0
end

"""
    effective_charge_distribution(model, masses) -> ChargeDistribution

The model reduced, at each heavy fragment mass, to the plain Gaussian of the same first two
moments: `ΔZ(A) = Z̄ − A Z_F/A_F` and the width `RMS(A)` on the lattice footing an evaluated
table carries. The even-odd structure survives as a ripple of period `≈ 2A_F/Z_F` in both. Meant
for comparison with tabulated distributions; the model itself serves the domain directly.
"""
function effective_charge_distribution(
    model::ZpModel,
    masses::AbstractVector{<:Integer};
    source::AbstractString = sprint(show, model),
)
    isempty(masses) && throw(ArgumentError("no masses to tabulate the model over"))
    ΔZ = Dict{Int, Float64}()
    σ_Z = Dict{Int, Float64}()
    for A in masses
        2 * A >= model.A_F ||
            throw(ArgumentError("tabulate heavy masses, A ≥ A_F/2; got $A"))
        Z̄, RMS = charge_moments(fragment_charge_yields(model, Int(A))...)
        ΔZ[Int(A)] = Z̄ - A * model.Z_F / model.A_F
        σ_Z[Int(A)] = RMS
    end
    return ChargeDistribution(
        ΔZ,
        σ_Z,
        DEFAULT_CHARGE_POLARIZATION,
        DEFAULT_CHARGE_DISPERSION,
        String(source),
    )
end

effective_charge_distribution(model::ZpModel, masses::AbstractRange{<:Integer}; kwargs...) =
    effective_charge_distribution(model, collect(masses); kwargs...)
