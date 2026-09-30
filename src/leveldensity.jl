"""
    LevelDensityModel

Systematics for the energy-independent level density parameter `a(A,Z)`. The sequential
equations close on the Fermi-gas relation `Ē_r = a T²`, which admits no explicit energy
dependence in `a`; phenomenological systematics are therefore used rather than a model with
`a(E*)`, such as the superfluid one, which would make the closure circular. Values are plain
`Float64`, since `a` is consumed by the TXE partition and the sequential solver, downstream of
the point where uncertainty propagation stops.
"""
abstract type LevelDensityModel end

"""
    LEVEL_DENSITY_MODELS

The level density models a configuration may select, in the spelling `[level_density] model`
takes: `"BSFG"` for [`BackShiftedFermiGas`](@ref) and `"GC"` for [`GilbertCameron`](@ref).
"""
const LEVEL_DENSITY_MODELS = ("BSFG", "GC")

"""
    ShellCorrectionTable

Gilbert-Cameron shell corrections `S(Z)` and `S(N)`, indexed by nucleon number.
"""
struct ShellCorrectionTable
    S_Z::Dict{Int, Float64}
    S_N::Dict{Int, Float64}
    source::String

    function ShellCorrectionTable(
        S_Z::Dict{Int, Float64},
        S_N::Dict{Int, Float64},
        source::AbstractString,
    )
        return new(S_Z, S_N, String(source))
    end
end

"""
    SHELL_CORRECTION_SPEC

Layout of the shell-correction file: one header row, then `n S_N S_Z`, the neutron correction
first, as in the file header `n S(N) S(Z)`. The Gilbert-Cameron expression evaluates `S(Z)` at
`n = Z` and `S(N)` at `n = N`, so the column order is not absorbed by the sum.
"""
const SHELL_CORRECTION_SPEC = TableSpec([:n, :S_N, :S_Z]; skip = 1)

"""
    read_shell_correction_table(path; spec = SHELL_CORRECTION_SPEC) -> ShellCorrectionTable

Read the Gilbert-Cameron shell corrections. Columns are taken by position, not by header text:
nucleon number, `S(N)`, `S(Z)`.
"""
function read_shell_correction_table(
    path::AbstractString;
    spec::TableSpec = SHELL_CORRECTION_SPEC,
)
    table = read_delimited_table(path, spec)
    nucleons = integer_column(table, :n)
    neutron = column(table, :S_N, Float64)
    proton = column(table, :S_Z, Float64)

    # Table III tabulates `S(Z)` to Z = 98 and `S(N)` to N = 150. The shipped file marks the
    # untabulated proton cells `NaN`; tables in circulation pad them with an exact 0.00 instead.
    # No tabulated `S(Z)` is zero, so trailing zeros are dropped as padding, and either way
    # `GilbertCameron` declines Z ≥ 99 rather than applying a correction the paper does not give.
    last_proton = something(findlast(value -> isfinite(value) && !iszero(value), proton), 0)

    protons = Dict{Int, Float64}()
    neutrons = Dict{Int, Float64}()
    for row in eachindex(nucleons)
        row <= last_proton &&
            isfinite(proton[row]) &&
            (protons[nucleons[row]] = proton[row])
        isfinite(neutron[row]) && (neutrons[nucleons[row]] = neutron[row])
    end
    return ShellCorrectionTable(protons, neutrons, String(path))
end

function Base.show(io::IO, corrections::ShellCorrectionTable)
    print(
        io,
        "ShellCorrectionTable(",
        length(corrections.S_Z),
        " entries from ",
        basename(corrections.source),
        ")",
    )
    return nothing
end

"""
    GilbertCameron(corrections)

Gilbert-Cameron systematics with shell corrections,

```
a(A,Z) = A [0.00917 (S(Z) + S(N)) + c],   c = 0.142 undeformed, 0.120 deformed
```

with `S(Z)` taken at the proton number and `S(N)` at the neutron number, both from Table III of
the paper. The offset follows eq. (20), `a/A = 0.00917S + 0.142`, for undeformed nuclei and
eq. (21), `a/A = 0.00917S + 0.120`, for deformed ones (p. 1457); [`is_deformed_gilbert_cameron`](@ref)
assigns the branch.

`deformed_branch = true`, the default, applies both branches as the paper prescribes; about a
third of the yield-weighted ²⁵²Cf fragment population lies in the first deformed region, where
the undeformed offset overestimates `a` by some 22 %. `deformed_branch = false` applies eq. (20)
everywhere, the convention of Tudora, Hambsch and Tobosaru, Eur. Phys. J. A 54, 87 (2018),
doi:10.1140/epja/i2018-12521-7, and of `FissionTemperatureRatio.jl`; it is required when an
`R_T(A_H)` produced by that package is inverted here with Gilbert-Cameron selected on both sides,
so that both form `a_L/a_H` the same way.

Source: A. Gilbert and A. G. W. Cameron, *Canadian Journal of Physics* **43**, 1446 (1965),
doi:10.1139/p65-139; eqs. (20) and (21) p. 1457, the deformation regions p. 1459, Table III
pp. 1453-1455.
"""
struct GilbertCameron <: LevelDensityModel
    corrections::ShellCorrectionTable
    deformed_branch::Bool

    function GilbertCameron(corrections::ShellCorrectionTable; deformed_branch::Bool = true)
        return new(corrections, deformed_branch)
    end
end

const GILBERT_CAMERON_SHELL_WEIGHT = 0.00917
const GILBERT_CAMERON_OFFSET_UNDEFORMED = 0.142
const GILBERT_CAMERON_OFFSET_DEFORMED = 0.120

"""
    is_deformed_gilbert_cameron(nuclide) -> Bool

Whether Gilbert-Cameron treats this nuclide as deformed, and so takes the `0.120` branch. The
regions are those of Gilbert and Cameron, Can. J. Phys. 43, 1446 (1965), doi:10.1139/p65-139,
p. 1459: "we shall consider the following to be regions of deformation: `54 ≤ Z ≤ 78,
86 ≤ N ≤ 122`" and "`86 ≤ Z ≤ 122, 130 ≤ N ≤ 182`", both bounds of a pair required together.
Light nuclei with `11 ≤ Z ≤ 19` follow the deformed line and `20 ≤ Z ≤ 29` the undeformed one
(p. 1464); these lie outside the fragment range. The first region covers the upper half of the
heavy fragment peak, `A ≳ 141`, in ²³⁵U(n,f).
"""
function is_deformed_gilbert_cameron(nuclide::Nuclide)
    Z = nuclide.Z
    N = neutron_number(nuclide)
    (54 <= Z <= 78 && 86 <= N <= 122) && return true
    (86 <= Z <= 122 && 130 <= N <= 182) && return true
    return 11 <= Z <= 19
end

"""
    BackShiftedFermiGas(table)

Von Egidy-Bucurescu systematics for the back-shifted Fermi gas,

```
a(A,Z) = (p₁ + p₂ δW) A^p₃
```

where the shell correction `δW = δW₀ + ½P_a′` is the difference `M_exp − M_LD` between the
experimental mass and a liquid-drop one, the order of eq. (7), plus a pairing term built from
the mass excesses of the neighbours `(A+2, Z+1)` and `(A−2, Z−1)`. Preferred over
Gilbert-Cameron for fission fragments: across the fragment range it stays close to the
superfluid model at fragment excitation energies, departing mainly near `A ≈ 130`, where the
`N = 82` and `Z = 50` shells close. The systematics is eq. (19) of the 2009 update, additive in
`S′`:

```
a = (p₁ + p₂ S′) A^p₃,   p₁ = 0.199(7), p₂ = 0.0096(4), p₃ = 0.869(7)
```

with the components

- `S(Z,N) = M_exp − M_LD` (2005, eq. 7), which in binding energies is `E_b^LD − E_b^exp`.
- `M_LD` from Pearson's Weizsäcker-type formula (2005, eq. 9), with `a_vol = −15.65`,
  `a_sf = 17.63`, `a_sym = 27.72`, `a_ss = −25.60` MeV and `r₀ = 1.233` fm, so the Coulomb
  coefficient is `3e²/5r₀ = 0.864/1.233` MeV with `e² = 1.44` MeV fm.
- `S′ = S + 0.5 P_a′` (2009, eq. 19) with
  `P_a′ = ½[M(A+2,Z+1) − 2M(A,Z) + M(A−2,Z−1)]` (2009, eq. 12), so the pairing enters as
  `¼[M(A+2,Z+1) − 2M(A,Z) + M(A−2,Z−1)]`.

`P_a′` carries no `(−1)^Z` alternation; it differs by `(−1)^(Z+1)` from the `P_a` of the 2005
paper and of Audi's tables. The 2009 form is used throughout, not the 2005 form
`a/A = p₁ + p₂ S′ + p₃ A`.

Only `a` is taken from the systematics. The papers fit it jointly with the back-shift
`E1 = −0.381 + 0.5 P_a′` (2009, eq. 20) for `ρ(U) ∝ exp(2√(a(U − E1)))`, whereas the sequential
equations close on the un-shifted `Ē_r = a T²`, as in eq. (1) of Tudora, Hambsch and Tobosaru,
Eur. Phys. J. A 54, 87 (2018), doi:10.1140/epja/i2018-12521-7, and in Eur. Phys. J. A 58, 126
(2022), doi:10.1140/epja/s10050-022-00766-y; `a` therefore enters as an effective parameter.
Over the fragment domain `E1` ranges from −1.7 to +1.9 MeV with a sign set by parity, so its
omission moves `T` by 1–3 % at the first
emission and by 18–40 % at the last sequence.

References: T. von Egidy and D. Bucurescu, *Phys. Rev. C* **72**, 044311 (2005),
doi:10.1103/PhysRevC.72.044311; erratum *Phys. Rev. C* **73**, 049901 (2006),
doi:10.1103/PhysRevC.73.049901; *Phys. Rev. C* **80**, 054310 (2009),
doi:10.1103/PhysRevC.80.054310.

"""
struct BackShiftedFermiGas <: LevelDensityModel
    mass_excess_table::MassExcessTable
    proton_excess::Float64
    neutron_excess::Float64

    function BackShiftedFermiGas(table::MassExcessTable)
        proton = mass_excess(table, Nuclide(1, 1))
        neutron = mass_excess(table, NEUTRON)
        (proton === nothing || neutron === nothing) && throw(
            ArgumentError(
                "the back-shifted Fermi gas systematics needs the mass excesses of ¹H and the " *
                "neutron, and $(table.source) supplies at least one of them not at all",
            ),
        )
        return new(table, value(proton), value(neutron))
    end
end

# Liquid-drop coefficients entering δW₀, in MeV.
const LDM_VOLUME = 15.65
const LDM_SURFACE = 17.63
const LDM_COULOMB = 0.864 / 1.233
const LDM_SYMMETRY_VOLUME = 27.72
const LDM_SYMMETRY_SURFACE = 25.6

# Coefficients of the level density parameter itself.
const BSFG_OFFSET = 1.99e-1
const BSFG_SHELL_WEIGHT = 9.6e-3
const BSFG_MASS_EXPONENT = 8.69e-1

"""
    liquid_drop_energy(model, nuclide) -> Float64

Liquid-drop binding energy in MeV, with a mass-dependent symmetry term.
"""
function liquid_drop_energy(::BackShiftedFermiGas, nuclide::Nuclide)
    A = nuclide.A
    Z = nuclide.Z
    symmetry = A * (LDM_SYMMETRY_VOLUME - LDM_SYMMETRY_SURFACE * A^(-1 / 3))
    η = (A - 2Z) / A
    return LDM_VOLUME * A - LDM_SURFACE * A^(2 / 3) - LDM_COULOMB * Z^2 * A^(-1 / 3) -
           symmetry * η^2
end

"""
    shell_correction(model, nuclide) -> Union{Float64,Nothing}

`S′ = S + ½P_a′` in MeV, the shell correction plus half the pairing term of the 2009
systematics,

```
½P_a′ = [Δ(A+2, Z+1) − 2Δ(A, Z) + Δ(A−2, Z−1)] / 4
```

`P_a′` is that of eq. (12) of von Egidy and Bucurescu, Phys. Rev. C 80, 054310 (2009),
doi:10.1103/PhysRevC.80.054310, consumed by its eq. (19); it carries no `(−1)^Z` alternation and
differs from the `P_d` of the 2005 paper and of Audi's tables. Returns `nothing` when a required
mass excess is unavailable, as for the most exotic fragments of a deterministic domain.
"""
function shell_correction(model::BackShiftedFermiGas, nuclide::Nuclide)
    Δ = mass_excess(model.mass_excess_table, nuclide)
    Δ === nothing && return nothing
    Δ_up = mass_excess(model.mass_excess_table, Nuclide(nuclide.Z + 1, nuclide.A + 2))
    Δ_up === nothing && return nothing
    nuclide.Z >= 1 && nuclide.A >= 2 || return nothing
    Δ_down = mass_excess(model.mass_excess_table, Nuclide(nuclide.Z - 1, nuclide.A - 2))
    Δ_down === nothing && return nothing

    experimental =
        nuclide.Z * model.proton_excess + neutron_number(nuclide) * model.neutron_excess -
        value(Δ)
    δW₀ = liquid_drop_energy(model, nuclide) - experimental
    pairing = (value(Δ_up) - 2 * value(Δ) + value(Δ_down)) / 4
    return δW₀ + pairing
end

"""
    level_density_parameter(model, nuclide) -> Union{Float64,Nothing}

The energy-independent level density parameter `a(A,Z)` in MeV⁻¹, or `nothing` when the
systematics cannot be evaluated or returns a non-positive value, which would give an imaginary
temperature through `Ē_r = a T²`.

# Examples

```jldoctest
julia> table = read_mass_excess_table(AME2020_MASS_EXCESS_FILE);

julia> a = level_density_parameter(BackShiftedFermiGas(table), Nuclide(52, 134));

julia> round(a; digits = 3)
7.707
```
"""
function level_density_parameter(model::BackShiftedFermiGas, nuclide::Nuclide)
    δW = shell_correction(model, nuclide)
    δW === nothing && return nothing
    a = (BSFG_OFFSET + BSFG_SHELL_WEIGHT * δW) * nuclide.A^BSFG_MASS_EXPONENT
    return a > 0 ? a : nothing
end

function level_density_parameter(model::GilbertCameron, nuclide::Nuclide)
    S_proton = get(model.corrections.S_Z, nuclide.Z, nothing)
    S_proton === nothing && return nothing
    S_neutron = get(model.corrections.S_N, neutron_number(nuclide), nothing)
    S_neutron === nothing && return nothing
    offset =
        model.deformed_branch && is_deformed_gilbert_cameron(nuclide) ?
        GILBERT_CAMERON_OFFSET_DEFORMED : GILBERT_CAMERON_OFFSET_UNDEFORMED
    a = nuclide.A * (GILBERT_CAMERON_SHELL_WEIGHT * (S_proton + S_neutron) + offset)
    return a > 0 ? a : nothing
end
