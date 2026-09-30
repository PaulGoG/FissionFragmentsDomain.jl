"""
    ChargeDistribution

Charge polarization `ΔZ(A)` and dispersion `σ_Z(A)` of the isobaric charge distribution, indexed
by fragment mass number, with the values used where a mass is not tabulated. The defaults
`⟨ΔZ⟩ = -0.5` and `⟨σ_Z⟩ = 0.6` are the conventional means applied when no `Zₚ` systematics is
available for the fission case. `ΔZ` is quoted for the heavy fragment and carries the opposite
sign for the light one.
"""
struct ChargeDistribution
    ΔZ::Dict{Int, Float64}
    σ_Z::Dict{Int, Float64}
    default_ΔZ::Float64
    default_σ_Z::Float64
    source::String

    function ChargeDistribution(
        ΔZ::Dict{Int, Float64},
        σ_Z::Dict{Int, Float64},
        default_ΔZ::Real,
        default_σ_Z::Real,
        source::AbstractString,
    )
        default_σ_Z > 0 || throw(
            ArgumentError(
                "the default charge dispersion must be positive, got $default_σ_Z",
            ),
        )
        return new(ΔZ, σ_Z, Float64(default_ΔZ), Float64(default_σ_Z), String(source))
    end
end

"""
    CHARGE_DISTRIBUTION_SPEC

Layout of the `ΔZ(A)`, `σ_Z(A)` files: one header row, then `A dZ sigma_Z`.
"""
const CHARGE_DISTRIBUTION_SPEC = TableSpec([:A, :dZ, :sigma_Z]; skip = 1)

const DEFAULT_CHARGE_POLARIZATION = -0.5
const DEFAULT_CHARGE_DISPERSION = 0.6

# `p(Z,A)` is the analytic Gaussian, not renormalized over the charges a domain retains, so a
# narrow charge window shows as lost weight. The yield construction `Y(A,Z,TKE) = p(Z,A) Y(A,TKE)`
# normalizes exactly, so no downstream quantity inherits the deviation of Σ_Z p from unity.

"""
    read_charge_distribution(path; spec, default_charge_polarization, default_charge_dispersion)
        -> ChargeDistribution

Read a tabulated `ΔZ(A)`, `σ_Z(A)`. Columns are taken by position, not by header text: mass
number, polarization, dispersion, so the reader is independent of the header text of files
written by other packages.
"""
function read_charge_distribution(
    path::AbstractString;
    spec::TableSpec = CHARGE_DISTRIBUTION_SPEC,
    default_charge_polarization::Real = DEFAULT_CHARGE_POLARIZATION,
    default_charge_dispersion::Real = DEFAULT_CHARGE_DISPERSION,
)
    table = read_delimited_table(path, spec)
    masses = integer_column(table, :A)
    polarization = column(table, :dZ, Float64)
    dispersion = column(table, :sigma_Z, Float64)

    ΔZ = Dict{Int, Float64}()
    σ_Z = Dict{Int, Float64}()
    for row in eachindex(masses)
        dispersion[row] > 0 || throw(
            ArgumentError(
                "$(table.path) row $row: the charge dispersion must be positive, got \
                 $(dispersion[row]) for A = $(masses[row])",
            ),
        )
        ΔZ[masses[row]] = polarization[row]
        σ_Z[masses[row]] = dispersion[row]
    end
    return ChargeDistribution(
        ΔZ,
        σ_Z,
        Float64(default_charge_polarization),
        Float64(default_charge_dispersion),
        String(path),
    )
end

"""
    mean_charge_distribution(; default_charge_polarization, default_charge_dispersion)
        -> ChargeDistribution

The fallback systematics with no tabulated values: every mass number takes the mean `ΔZ` and
`σ_Z`.
"""
function mean_charge_distribution(;
    default_charge_polarization::Real = DEFAULT_CHARGE_POLARIZATION,
    default_charge_dispersion::Real = DEFAULT_CHARGE_DISPERSION,
)
    return ChargeDistribution(
        Dict{Int, Float64}(),
        Dict{Int, Float64}(),
        Float64(default_charge_polarization),
        Float64(default_charge_dispersion),
        "mean values",
    )
end

function Base.show(io::IO, distribution::ChargeDistribution)
    print(
        io,
        "ChargeDistribution(",
        length(distribution.ΔZ),
        " masses from ",
        basename(distribution.source),
        ")",
    )
    return nothing
end

"""
    charge_polarization(distribution, A) -> Float64

Charge polarization `ΔZ` for heavy-fragment mass `A`, falling back to the mean value.
"""
charge_polarization(distribution::ChargeDistribution, A::Integer) =
    get(distribution.ΔZ, Int(A), distribution.default_ΔZ)

"""
    charge_dispersion(distribution, A) -> Float64

Dispersion `σ_Z` of the isobaric charge distribution at mass `A`, falling back to the mean value.
"""
charge_dispersion(distribution::ChargeDistribution, A::Integer) =
    get(distribution.σ_Z, Int(A), distribution.default_σ_Z)

"""
    unchanged_charge_distribution(system, A) -> Float64

`Z_UCD(A) = A Z₀ / A₀`, the charge a fragment would carry if the compound nucleus split
without charge rearrangement.
"""
function unchanged_charge_distribution(system::FissioningSystem, A::Integer)
    return A * system.compound.Z / system.compound.A
end

"""
    most_probable_charge(system, A, ΔZ) -> Float64

`Zₚ(A) = Z_UCD(A) + ΔZ(A)`, the centre of the isobaric charge distribution.
"""
function most_probable_charge(system::FissioningSystem, A::Integer, ΔZ::Real)
    return unchanged_charge_distribution(system, A) + ΔZ
end

"""
    charge_probability(Z, Zₚ, σ_Z) -> Float64

The isobaric charge distribution, a Gaussian of dispersion `σ_Z` centred on `Zₚ`:

```
p(Z, A) = exp(−(Z − Zₚ)² / 2 σ_Z²) / (√(2π) σ_Z)
```
"""
function charge_probability(Z::Real, Zₚ::Real, σ_Z::Real)
    σ_Z > 0 || throw(ArgumentError("the charge dispersion must be positive, got $σ_Z"))
    return exp(-(Z - Zₚ)^2 / (2 * σ_Z^2)) / (sqrt(2π) * σ_Z)
end

"""
    charge_numbers(Zₚ, count) -> UnitRange{Int}

The `count` integer charge numbers nearest `Zₚ`, centred on it. `count` must be odd, so that
the set is symmetric about the rounded most probable charge. Rounding follows Julia's default,
ties to even: a `Zₚ` exactly halfway between two integers, as at symmetric fission under the
mean polarization (`Zₚ = Z₀/2 − 0.5`), rounds to the even neighbour.
"""
function charge_numbers(Zₚ::Real, count::Integer)
    isodd(count) ||
        throw(ArgumentError("the number of charges per mass must be odd, got $count"))
    count > 0 || throw(ArgumentError("the number of charges per mass must be positive"))
    centre = round(Int, Zₚ)
    half = (count - 1) ÷ 2
    return (centre - half):(centre + half)
end
