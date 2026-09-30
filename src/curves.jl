"""
    SegmentedCurve(masses, values, source = "")

A piecewise-linear function of heavy-fragment mass, given as segment endpoints `(A_H, y)` in
ascending `A_H`. `E*_H/TXE (A_H)` and `R_T (A_H)` are supplied in this form. `source` records
the file the curve was read from and enters the run metadata. Outside the tabulated range the
function evaluates to `nothing` on both sides, so a fit is never extrapolated.
"""
struct SegmentedCurve
    masses::Vector{Int}
    values::Vector{Float64}
    source::String

    function SegmentedCurve(
        masses::Vector{Int},
        values::Vector{Float64},
        source::AbstractString = "",
    )
        length(masses) == length(values) ||
            throw(ArgumentError("segment masses and values must have equal length"))
        length(masses) >= 2 ||
            throw(ArgumentError("a segmented curve needs at least two points"))
        issorted(masses) || throw(
            ArgumentError("segment points must be given in ascending mass, got $masses"),
        )
        allunique(masses) ||
            throw(ArgumentError("segment points must have distinct masses, got $masses"))
        return new(copy(masses), copy(values), String(source))
    end
end

"""
    SegmentedCurve(points::AbstractVector{<:Tuple})

Build from `(A_H, value)` pairs.

# Examples

```jldoctest
julia> curve = SegmentedCurve([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]);

julia> curve(130)
0.14

julia> round(curve(140); digits = 4)
0.3735
```
"""
function SegmentedCurve(
    points::AbstractVector{<:Tuple{Integer, Real}},
    source::AbstractString = "",
)
    return SegmentedCurve(
        [Int(p[1]) for p in points],
        [Float64(p[2]) for p in points],
        source,
    )
end

Base.length(curve::SegmentedCurve) = length(curve.masses)

function Base.show(io::IO, curve::SegmentedCurve)
    print(
        io,
        "SegmentedCurve(",
        length(curve),
        " points, A_H ∈ ",
        first(curve.masses),
        ":",
        last(curve.masses),
        ")",
    )
    return nothing
end

"""
    (curve)(A) -> Union{Float64,Nothing}

Evaluate the curve at heavy-fragment mass `A`, or `nothing` outside the tabulated range or where
the interpolated value is non-positive.
"""
function (curve::SegmentedCurve)(A::Integer)
    A < first(curve.masses) && return nothing
    for index in 2:length(curve)
        if curve.masses[index] >= A
            x₀ = curve.masses[index - 1]
            x₁ = curve.masses[index]
            y₀ = curve.values[index - 1]
            y₁ = curve.values[index]
            slope = (y₁ - y₀) / (x₁ - x₀)
            y = y₀ + slope * (A - x₀)
            return y > 0 ? y : nothing
        end
    end
    return nothing
end

"""
    TEMPERATURE_RATIO_CURVE_SPEC

Layout of a tabulated `R_T(A_H)`: one header row, then `A_H,R_T` and an optional
`R_T_uncertainty`, comma-delimited, as written by a temperature-ratio run. The uncertainty is
read and discarded: the curve feeds an excitation-energy partition that carries no uncertainty
downstream.
"""
const TEMPERATURE_RATIO_CURVE_SPEC = TableSpec(
    [:A_H, :R_T, :R_T_uncertainty];
    delimiter = ',',
    skip = 1,
    optional = (:R_T_uncertainty,),
)

"""
    EXCITATION_RATIO_CURVE_SPEC

Layout of a tabulated `E*_H/TXE (A_H)`: one header row, then `A_H,excitation_ratio` and an
optional `excitation_ratio_uncertainty`, comma-delimited. Counterpart of
[`TEMPERATURE_RATIO_CURVE_SPEC`](@ref) for the excitation-energy ratio.
"""
const EXCITATION_RATIO_CURVE_SPEC = TableSpec(
    [:A_H, :excitation_ratio, :excitation_ratio_uncertainty];
    delimiter = ',',
    skip = 1,
    optional = (:excitation_ratio_uncertainty,),
)

"""
    read_segmented_curve(path; spec = TEMPERATURE_RATIO_CURVE_SPEC) -> SegmentedCurve

Read a curve tabulated against heavy-fragment mass. Columns are taken by position, not by
header text: the abscissa is the first column of `spec` and the ordinate the second.

The format is that of a temperature-ratio extraction, which tabulates `R_T = T_L/T_H` at every
integer `A_H` of its fit range, so [`SegmentedCurve`](@ref) evaluates only tabulated points and
no interpolation occurs. The breakpoints of the extraction's piecewise-linear fit of `r_ν` are
not equivalent, because `R_T = √[(1−r_ν)/(R_a r_ν)]` inherits the shell structure of
`R_a(A_H)`; see [`temperature_ratio_path`](@ref).

Each curve spans its own fit range and evaluates to `nothing` outside it. A partition that
applies per-`(A_H, Z_H)` level density parameters and reduces over `Z` afterwards corresponds to
`⟨a_L/a_H⟩` rather than `⟨a_L⟩/⟨a_H⟩`, so its `R_T` should come from a run with
`ratio_averaging = "mean_of_ratios"`; the setting is recorded in the run identifier the file name
carries.
"""
function read_segmented_curve(
    path::AbstractString;
    spec::TableSpec = TEMPERATURE_RATIO_CURVE_SPEC,
)
    table = read_delimited_table(path, spec)
    masses = integer_column(table, first(spec.columns))
    values = column(table, spec.columns[2], Float64)
    return SegmentedCurve(masses, values, path)
end
