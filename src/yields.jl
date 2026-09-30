"""
    MASS_ENERGY_YIELD_SPEC

Layout of the measured pre-neutron mass yields: one header row, then `A TKE Y` and an optional
`Y_uncertainty`.
"""
const MASS_ENERGY_YIELD_SPEC =
    TableSpec([:A, :TKE, :Y, :Y_uncertainty]; skip = 1, optional = (:Y_uncertainty,))

"""
    MassEnergyYield

Measured pre-neutron fragment mass yield `Y(A, TKE)`, on the grid the measurement was reported
on, which may differ from the model's `TKE` grid; [`coverage`](@ref) quantifies the mismatch.
`A` runs over both fragments, light and heavy, so the distribution sums to 200 % rather than
100 %.
"""
struct MassEnergyYield
    values::Dict{Tuple{Int, Float64}, Float64}
    masses::Vector{Int}
    energies::Vector{Float64}
    source::String
end

Base.length(Y::MassEnergyYield) = length(Y.values)

function Base.show(io::IO, Y::MassEnergyYield)
    print(
        io,
        "MassEnergyYield(",
        length(Y.masses),
        " masses × ",
        length(Y.energies),
        " TKE, ",
        length(Y.values),
        " cells, from ",
        basename(Y.source),
        ")",
    )
    return nothing
end

"""
    read_mass_energy_yield(path; spec = MASS_ENERGY_YIELD_SPEC) -> MassEnergyYield

Read a measured `Y(A, TKE)` distribution. Columns are taken by position, not by header text:
mass, total kinetic energy, yield, uncertainty. Cells with a non-finite yield are dropped, since
an unmeasured `(A, TKE)` is marked `NaN` and is not a yield of zero. A repeated `(A, TKE)` with
identical yields, such as `A = A₀/2` listed from each side of the split, is counted once; a
repeat with differing yields throws an `ArgumentError`.
"""
function read_mass_energy_yield(
    path::AbstractString;
    spec::TableSpec = MASS_ENERGY_YIELD_SPEC,
)
    table = read_delimited_table(path, spec)
    masses = integer_column(table, :A)
    energies = column(table, :TKE, Float64)
    measured = column(table, :Y, Float64)

    values = Dict{Tuple{Int, Float64}, Float64}()
    for row in eachindex(masses)
        isfinite(measured[row]) || continue
        key = (masses[row], energies[row])
        if haskey(values, key)
            values[key] == measured[row] || throw(
                ArgumentError(
                    "$(table.path) row $row: $(key) appears more than once with different " *
                    "yields, $(values[key]) and $(measured[row])",
                ),
            )
            continue
        end
        values[key] = measured[row]
    end
    isempty(values) && throw(ArgumentError("$(table.path) holds no finite yields"))
    return MassEnergyYield(
        values,
        sort(unique(masses)),
        sort(unique(energies)),
        String(path),
    )
end

"""
    MASS_YIELD_SPEC
    MEAN_KINETIC_ENERGY_SPEC

Layouts of the two marginals a factorized yield is built from:
`A Y [Y_uncertainty]` and `A TKE [TKE_uncertainty]`, each with one header row, as written by a
retrieval for the `Y_vs_A` and `TKE_vs_A` observables.
"""
const MASS_YIELD_SPEC =
    TableSpec([:A, :Y, :Y_uncertainty]; skip = 1, optional = (:Y_uncertainty,))

@doc (@doc MASS_YIELD_SPEC)
const MEAN_KINETIC_ENERGY_SPEC =
    TableSpec([:A, :TKE, :TKE_uncertainty]; skip = 1, optional = (:TKE_uncertainty,))

"""
    KINETIC_ENERGY_DISPERSION_SPEC

Layout of a tabulated `σ_TKE(A)`: one header row, then `A sigma_TKE`, in MeV.
"""
const KINETIC_ENERGY_DISPERSION_SPEC = TableSpec([:A, :sigma_TKE]; skip = 1)

"""
    KineticEnergyDispersion

The dispersion `σ_TKE(A)` of `P(TKE|A)` in MeV, indexed by fragment mass, with `default`, a
single width, used where a mass is not tabulated. Structured as [`ChargeDistribution`](@ref).
"""
struct KineticEnergyDispersion
    σ_TKE::Dict{Int, Float64}
    default::Float64
    source::String

    function KineticEnergyDispersion(
        σ_TKE::Dict{Int, Float64},
        default::Real,
        source::AbstractString,
    )
        default > 0 || throw(
            ArgumentError(
                "the default kinetic energy dispersion must be positive, got $default",
            ),
        )
        return new(σ_TKE, Float64(default), String(source))
    end
end

"""
    DEFAULT_KINETIC_ENERGY_DISPERSION

The single width of `P(TKE|A)` used where no table applies, in MeV. A stated value, not a fit;
the yield-weighted widths of the shipped joint distributions lie between 7.9 and 9.6 MeV.
"""
const DEFAULT_KINETIC_ENERGY_DISPERSION = 10.0

"""
    read_kinetic_energy_dispersion(path; spec, default) -> KineticEnergyDispersion

Read a tabulated `σ_TKE(A)`. Columns are taken by position, not by header text: mass number,
then dispersion.
"""
function read_kinetic_energy_dispersion(
    path::AbstractString;
    spec::TableSpec = KINETIC_ENERGY_DISPERSION_SPEC,
    default::Real = DEFAULT_KINETIC_ENERGY_DISPERSION,
)
    table = read_delimited_table(path, spec)
    masses = integer_column(table, :A)
    dispersions = column(table, :sigma_TKE, Float64)
    σ_TKE = Dict{Int, Float64}()
    for row in eachindex(masses)
        dispersions[row] > 0 || throw(
            ArgumentError(
                "$(table.path) row $row: the kinetic energy dispersion must be positive, got \
                 $(dispersions[row]) for A = $(masses[row])",
            ),
        )
        σ_TKE[masses[row]] = dispersions[row]
    end
    return KineticEnergyDispersion(σ_TKE, default, String(path))
end

"""
    uniform_kinetic_energy_dispersion(σ_TKE) -> KineticEnergyDispersion

One width at every mass, the fallback of a reconstruction without a table.
"""
uniform_kinetic_energy_dispersion(σ_TKE::Real = DEFAULT_KINETIC_ENERGY_DISPERSION) =
    KineticEnergyDispersion(Dict{Int, Float64}(), σ_TKE, "uniform")

"""
    kinetic_energy_dispersion(dispersion, A) -> Float64

`σ_TKE` at fragment mass `A`, falling back to the single value where the table does not reach.
"""
kinetic_energy_dispersion(dispersion::KineticEnergyDispersion, A::Integer) =
    get(dispersion.σ_TKE, Int(A), dispersion.default)

function Base.show(io::IO, dispersion::KineticEnergyDispersion)
    print(
        io,
        "KineticEnergyDispersion(",
        length(dispersion.σ_TKE),
        " masses, default ",
        dispersion.default,
        " MeV, from ",
        basename(dispersion.source),
        ")",
    )
    return nothing
end

"""
    factorized_yield(mass_yield_path, kinetic_energy_path, energies; dispersion) -> MassEnergyYield

Reconstruct `Y(A, TKE)` from the two marginals that were measured, where the joint distribution
was not:

```
Y(A, TKE) = Y(A) · G(TKE; ⟨TKE⟩(A), σ)
```

with `G` a Gaussian normalized over `energies`, from a measured pre-neutron mass yield `Y(A)` and
a measured mean total kinetic energy `⟨TKE⟩(A)`. The result is a reconstruction that retains the
mass dependence of `⟨TKE⟩`; a single mean at every mass over-subtracts at the asymmetric
splits, where `⟨TKE⟩` is about 20 MeV lower than in the peaks, and inverts the heavy branch of
`ν(A)`. The mass correlation of the width is carried only as far as `dispersion` supplies it,
either a measured [`KineticEnergyDispersion`](@ref) or a single width in MeV.

A mass without a tabulated `⟨TKE⟩` is dropped rather than extrapolated. A complementary mass,
`compound_mass − A`, not tabulated itself inherits the mean of its partner, since `TKE` is a
property of the split.

With `symmetrize = true` the exact pre-neutron identities `Y(A) = Y(A₀ − A)`,
`⟨TKE⟩(A) = ⟨TKE⟩(A₀ − A)` and `σ_TKE(A) = σ_TKE(A₀ − A)` are imposed: the two fragments of a
split are counted in one event, so their masses have one yield and one kinetic energy.
- Where both complements are tabulated, each takes their mean. Measured values that differ, such
  as a backing loss on one side of a double-energy measurement, then enter once and alike.
- Where only one is tabulated, it serves the other as before.
Pre-neutron data should be symmetrized. The default, `false`, keeps the tables as measured.
The result is renormalized to `total`. Columns are taken by position in both files, not by
header text: mass and yield in the first, mass and mean kinetic energy in the second.
"""
function factorized_yield(
    mass_yield_path::AbstractString,
    kinetic_energy_path::AbstractString,
    energies,
    compound_mass::Integer;
    dispersion::Union{Real, KineticEnergyDispersion} = DEFAULT_KINETIC_ENERGY_DISPERSION,
    total::Real = DEFAULT_YIELD_TOTAL,
    symmetrize::Bool = false,
)
    widths =
        dispersion isa KineticEnergyDispersion ? dispersion :
        uniform_kinetic_energy_dispersion(dispersion)
    grid = sort(unique(Float64.(collect(energies))))
    isempty(grid) && throw(ArgumentError("the kinetic energy grid is empty"))

    mass_table = read_delimited_table(mass_yield_path, MASS_YIELD_SPEC)
    masses = integer_column(mass_table, :A)
    yields = column(mass_table, :Y, Float64)

    energy_table = read_delimited_table(kinetic_energy_path, MEAN_KINETIC_ENERGY_SPEC)
    energy_masses = integer_column(energy_table, :A)
    means = column(energy_table, :TKE, Float64)

    # TKE belongs to the split, so a mean tabulated against the heavy mass serves the light
    # partner too.
    measured_mean = Dict{Int, Float64}()
    for row in eachindex(energy_masses)
        isfinite(means[row]) && means[row] > 0 || continue
        measured_mean[energy_masses[row]] = means[row]
    end
    symmetrize && (measured_mean = _complement_averaged(measured_mean, compound_mass))
    mean_by_mass = copy(measured_mean)
    for (A, mean) in measured_mean
        complement = Int(compound_mass) - A
        haskey(mean_by_mass, complement) || (mean_by_mass[complement] = mean)
    end
    isempty(mean_by_mass) &&
        throw(ArgumentError("$(energy_table.path) holds no positive mean kinetic energies"))

    measured_yield = Dict{Int, Float64}()
    for row in eachindex(masses)
        isfinite(yields[row]) &&
            yields[row] > 0 &&
            (measured_yield[masses[row]] = yields[row])
    end
    symmetrize && (measured_yield = _complement_averaged(measured_yield, compound_mass))
    width_by_mass =
        symmetrize ? _complement_averaged(widths.σ_TKE, compound_mass) : widths.σ_TKE

    values = Dict{Tuple{Int, Float64}, Float64}()
    kept = Int[]
    running = 0.0
    for A in sort(collect(keys(measured_yield)))
        Y = measured_yield[A]
        mean = get(mean_by_mass, A, nothing)
        mean === nothing && continue
        σ_TKE =
            symmetrize ? get(width_by_mass, A) do
                get(width_by_mass, Int(compound_mass) - A, widths.default)
            end : kinetic_energy_dispersion(widths, A)
        shape = [exp(-(TKE - mean)^2 / (2 * σ_TKE^2)) for TKE in grid]
        norm = sum(shape)
        norm > 0 || continue
        push!(kept, A)
        for (index, TKE) in pairs(grid)
            cell = Y * shape[index] / norm
            values[(A, TKE)] = cell
            running += cell
        end
    end
    isempty(values) &&
        throw(ArgumentError("no mass carries both a yield and a mean kinetic energy; \
                       $(mass_table.path) and $(energy_table.path) do not overlap"))

    # Renormalized to the 200 % of a measured two-fragment distribution, so factorized and
    # measured yields weight a run alike.
    scale = total / running
    for key in keys(values)
        values[key] *= scale
    end
    return MassEnergyYield(
        values,
        sort(kept),
        grid,
        string(basename(mass_yield_path), " × ", basename(kinetic_energy_path)),
    )
end

# The mean of the values at A and at its complement where both are tabulated; a value whose
# complement is absent is kept as it is.
function _complement_averaged(values::AbstractDict{Int, Float64}, compound_mass::Integer)
    averaged = Dict{Int, Float64}()
    for (A, value) in values
        partner = get(values, Int(compound_mass) - A, nothing)
        averaged[A] = partner === nothing ? value : (value + partner) / 2
    end
    return averaged
end

"""
    symmetrized_yield(yield, compound_mass) -> MassEnergyYield

`yield` with the exact pre-neutron identity `Y(A, TKE) = Y(A₀ − A, TKE)` imposed, `A₀` being
`compound_mass`. The two fragments of a split share one event, and so one yield at each total
kinetic energy.
- Where both cells `(A, TKE)` and `(A₀ − A, TKE)` are measured, each takes their mean.
- A cell whose complement is not measured is kept as it is.
The normalization is unchanged: the mean of two cells preserves their sum.
"""
function symmetrized_yield(yield::MassEnergyYield, compound_mass::Integer)
    values = Dict{Tuple{Int, Float64}, Float64}()
    for ((A, TKE), Y) in yield.values
        partner = get(yield.values, (Int(compound_mass) - A, TKE), nothing)
        values[(A, TKE)] = partner === nothing ? Y : (Y + partner) / 2
    end
    return MassEnergyYield(values, yield.masses, yield.energies, yield.source)
end

"""
    fragment_charge_probabilities(domain, distribution) -> Dict{Nuclide,Float64}

The isobaric charge probability `p(Z,A)` of every fragment the domain produces, the factor of
`Y(A,Z,TKE) = p(Z,A) Y(A,TKE)`. Each nuclide receives one value, evaluated against the `Zₚ` of
its own mass, neither accumulated over the fragmentations it appears in nor inherited from a
complementary fragment; the distinction matters only at `A = A₀/2`. This differs from
[`fragments`](@ref), which sums a nuclide's contributions over the domain. `p` is not
renormalized over the retained charges; see [`fragment_charge_probability`](@ref).
"""
function fragment_charge_probabilities(domain::FragmentationDomain, model::ChargeModel)
    probabilities = Dict{Nuclide, Float64}()
    for entry in domain.entries
        for fragment in (entry.light, entry.heavy)
            haskey(probabilities, fragment) && continue
            probabilities[fragment] =
                fragment_charge_probability(model, domain.system, fragment)
        end
    end
    return probabilities
end

"""
    FragmentYield

The yield of every fragment configuration, `Y(A, Z, TKE)`, normalized so the whole distribution
sums to `total`.
"""
struct FragmentYield
    values::Dict{Tuple{Int, Int, Float64}, Float64}
    total::Float64
end

Base.length(Y::FragmentYield) = length(Y.values)
Base.sum(Y::FragmentYield) = sum(values(Y.values))

function Base.show(io::IO, Y::FragmentYield)
    print(
        io,
        "FragmentYield(",
        length(Y.values),
        " cells, Σ = ",
        round(Y.total; digits = 4),
        ")",
    )
    return nothing
end

"""
    DEFAULT_YIELD_TOTAL

The normalization of a fragment yield distribution: 200, since a fission event produces two
fragments and each mass distribution is quoted as a percentage.
"""
const DEFAULT_YIELD_TOTAL = 200.0

"""
    build_fragment_yield(domain, distribution, experimental; total = DEFAULT_YIELD_TOTAL)
        -> FragmentYield

Spread the measured `Y(A, TKE)` over charge numbers,

```
Y(A, Z, TKE) = p(Z, A) Y(A, TKE)
```

and renormalize the result to `total`. The weight not captured by the retained charge numbers
is divided out once over the whole distribution, not per mass row, so the unnormalized
`p(Z,A)` does not distort the shape. Cells are taken on the measurement's `TKE` grid, not the
model's; [`coverage`](@ref) quantifies the overlap.
"""
function build_fragment_yield(
    domain::FragmentationDomain,
    distribution::ChargeModel,
    experimental::MassEnergyYield;
    total::Real = DEFAULT_YIELD_TOTAL,
)
    total > 0 ||
        throw(ArgumentError("the yield normalization must be positive, got $total"))
    probabilities = fragment_charge_probabilities(domain, distribution)

    values = Dict{Tuple{Int, Int, Float64}, Float64}()
    running = 0.0
    for (nuclide, p) in probabilities, TKE in experimental.energies
        measured = get(experimental.values, (nuclide.A, TKE), nothing)
        measured === nothing && continue
        contribution = p * measured
        values[(nuclide.A, nuclide.Z, TKE)] = contribution
        running += contribution
    end
    isempty(values) &&
        throw(ArgumentError("no fragment of the domain has a measured yield"))
    running > 0 || throw(ArgumentError("the fragment yields sum to $running"))

    scale = total / running
    for key in keys(values)
        values[key] *= scale
    end
    return FragmentYield(values, Float64(total))
end

"""
    coverage(yields, energies) -> NamedTuple

The part of a yield distribution a model sweep reaches: `reached`, the fraction of the total
yield on `energies`; `missed`, the remainder; `swept` and `tabulated`, the two `TKE` grids. A
sweep coarser than the measurement averages over part of the distribution only.

# Examples

```jldoctest
julia> using FissionFragmentsDomain

julia> yields = FragmentYield(Dict((140, 54, 170.0) => 3.0, (140, 54, 171.0) => 1.0), 4.0);

julia> report = coverage(yields, 170.0:2.0:174.0);

julia> report.reached
0.75
```
"""
function coverage(yields::FragmentYield, energies)
    wanted = Set(Float64.(energies))
    on = 0.0
    off = 0.0
    present = Set{Float64}()
    for ((_, _, TKE), Y) in yields.values
        push!(present, TKE)
        TKE in wanted ? (on += Y) : (off += Y)
    end
    whole = on + off
    return (
        reached = whole == 0 ? 0.0 : on / whole,
        missed = whole == 0 ? 0.0 : off / whole,
        swept = sort(collect(wanted)),
        tabulated = sort(collect(present)),
    )
end

"""
    yield_over(yields, dimension) -> Vector{Tuple{Float64,Float64}}

Project the yield onto one of its arguments, `:mass`, `:charge` or `:kinetic_energy`, summing
over the other two. Sorted by the argument.

Each projection sums to the same total as the distribution it came from.
"""
function yield_over(yields::FragmentYield, dimension::Symbol)
    index = if dimension === :mass
        1
    elseif dimension === :charge
        2
    elseif dimension === :kinetic_energy
        3
    else
        throw(
            ArgumentError(
                "unknown yield dimension $dimension; expected :mass, :charge or :kinetic_energy",
            ),
        )
    end
    projected = Dict{Float64, Float64}()
    for (key, Y) in yields.values
        argument = Float64(key[index])
        projected[argument] = get(projected, argument, 0.0) + Y
    end
    return sort!([(argument, Y) for (argument, Y) in projected]; by = first)
end

"""
    mean_kinetic_energy_by_mass(yields) -> Vector{Tuple{Int,Float64}}

`⟨TKE⟩(A)`, the yield-weighted mean total kinetic energy at each fragment mass.

Defined at every mass the distribution covers, light and heavy; the two halves carry the same
information, since `TKE` belongs to the pair.
"""
function mean_kinetic_energy_by_mass(yields::FragmentYield)
    weighted = Dict{Int, Float64}()
    total = Dict{Int, Float64}()
    for ((A, _, TKE), Y) in yields.values
        weighted[A] = get(weighted, A, 0.0) + Y * TKE
        total[A] = get(total, A, 0.0) + Y
    end
    return sort!(
        [(A, weighted[A] / total[A]) for A in keys(total) if total[A] > 0];
        by = first,
    )
end

"""
    MassYield

A pre-neutron fragment mass yield `Y(A)` with its uncertainties, ascending in `A`, as a
temperature-ratio extraction weights its total average. `σY` holds `missing` where the source
quotes no uncertainty, never zero. No normalisation is imposed: an average divides by the weights
it used, so a distribution normalised to one peak, to 2, or to 200 % gives the same answer.
`label` identifies the distribution in tables and legends; `source` is the file it was read from.
"""
struct MassYield
    A::Vector{Int}
    Y::Vector{Float64}
    σY::Vector{Union{Missing, Float64}}
    label::String
    source::String
end

Base.length(yields::MassYield) = length(yields.A)

"""
    mass_yield(yields, A) -> Union{Tuple{Float64,Union{Missing,Float64}},Nothing}

Yield and uncertainty at mass `A`, or `nothing` where the distribution has no entry.
"""
function mass_yield(yields::MassYield, A::Integer)
    index = findfirst(==(A), yields.A)
    return index === nothing ? nothing : (yields.Y[index], yields.σY[index])
end

"""
    read_mass_yield(path; label = "", spec = MASS_YIELD_SPEC) -> MassYield

Read `Y(A)` laid out as `A Y [Y_uncertainty]` with one header row, columns taken by position. A
non-finite or non-positive uncertainty is stored as `missing`. A negative yield or a repeated mass
throws an `ArgumentError` naming the file. `label` defaults to the file name without extension.

# Examples

```jldoctest
julia> path = tempname();

julia> write(path, "A Y Y_uncertainty\\n140 6.1 0.2\\n132 4.8 0\\n");

julia> yields = read_mass_yield(path; label = "example");

julia> yields.A, yields.Y
([132, 140], [4.8, 6.1])

julia> mass_yield(yields, 132)
(4.8, missing)
```
"""
function read_mass_yield(
    path::AbstractString;
    label::AbstractString = "",
    spec::TableSpec = MASS_YIELD_SPEC,
)
    table = read_delimited_table(path, spec)
    masses = integer_column(table, spec.columns[1])
    measured = column(table, spec.columns[2], Float64)
    quoted = column_or(table, spec.columns[3], Float64, NaN)

    allunique(masses) || throw(ArgumentError("$(table.path) repeats a mass number"))
    for row in eachindex(masses)
        measured[row] >= 0 ||
            throw(ArgumentError("$(table.path) has a negative yield at A = $(masses[row])"))
    end
    σY = Union{Missing, Float64}[isfinite(σ) && σ > 0 ? σ : missing for σ in quoted]
    order = sortperm(masses)
    name = isempty(label) ? first(splitext(basename(path))) : String(label)
    return MassYield(masses[order], measured[order], σY[order], name, String(path))
end
