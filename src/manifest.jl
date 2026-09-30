"""
    ManifestCurve

One `[[segmented_curve]]` of a temperature-ratio manifest: its label, its kind, and the two files
the producing run tabulated it in. It describes a curve; the function itself is a
[`SegmentedCurve`](@ref) read from one of these files. `kind` is `"dataset"` for a curve
extracted from a single measurement and `"systematic_trend"` for the combined trend fitted
across them. Only `temperature_ratio_file` is read, through [`temperature_ratio_path`](@ref);
`multiplicity_ratio_pivots_file` holds `r_ν`, not `R_T`, and is recorded only.
"""
struct ManifestCurve
    label::String
    kind::String
    temperature_ratio_file::String
    multiplicity_ratio_pivots_file::String

    function ManifestCurve(
        label::AbstractString,
        kind::AbstractString,
        temperature_ratio_file::AbstractString,
        multiplicity_ratio_pivots_file::AbstractString,
    )
        return new(
            String(label),
            String(kind),
            String(temperature_ratio_file),
            String(multiplicity_ratio_pivots_file),
        )
    end
end

function Base.show(io::IO, curve::ManifestCurve)
    print(io, "ManifestCurve(", repr(curve.label), ", ", curve.kind, ")")
    return nothing
end

"""
    MANIFEST_CURVE_KINDS

What a `[[segmented_curve]]` may declare itself to be: extracted from one dataset, or the
systematic trend combined across them.
"""
const MANIFEST_CURVE_KINDS = ("dataset", "systematic_trend")

"""
    SYSTEMATIC_TREND_LABEL

The label of the systematic-trend curve in a temperature-ratio manifest, `"systematic_trend"`,
the same token as its kind. A label is an identifier that a consumer's configuration names and
that file names carry, so it contains no whitespace; producers write this constant rather than a
spelling of their own.
"""
const SYSTEMATIC_TREND_LABEL = "systematic_trend"

"""
    MANIFEST_ORDINATE
    MANIFEST_ABSCISSA

The ordinate and abscissa a temperature-ratio manifest must declare: `R_T` against `A_H`. A
manifest declaring any other quantity is refused, since `r_ν` and similar quantities occupy the
same column position and would otherwise be read without error.
"""
const MANIFEST_ORDINATE = "R_T"

@doc (@doc MANIFEST_ORDINATE)
const MANIFEST_ABSCISSA = ["A_H"]

"""
    ManifestSystem

The `[system]` table of a manifest: the fissioning system the producing run extracted its curves
from. Field names mirror the keys. It is recorded, not matched against the consuming run's
system, since a run may take its curve from another system's manifest, e.g. a resonance-channel
run using the thermal system's trend.
"""
struct ManifestSystem
    label::String
    notation::String
    target_A::Int
    target_Z::Int
    channel::String
    reaction::String
    incident_energy_MeV::Float64
    compound_A::Int
    compound_Z::Int

    # Explicit inner constructor: concrete conversions at the call site, no generic `convert`
    # path into the type.
    function ManifestSystem(
        label::AbstractString,
        notation::AbstractString,
        target_A::Integer,
        target_Z::Integer,
        channel::AbstractString,
        reaction::AbstractString,
        incident_energy_MeV::Real,
        compound_A::Integer,
        compound_Z::Integer,
    )
        return new(
            String(label),
            String(notation),
            Int(target_A),
            Int(target_Z),
            String(channel),
            String(reaction),
            Float64(incident_energy_MeV),
            Int(compound_A),
            Int(compound_Z),
        )
    end
end

function Base.show(io::IO, system::ManifestSystem)
    print(io, "ManifestSystem(", system.label, ")")
    return nothing
end

"""
    TemperatureRatioManifest

The run record a temperature-ratio extraction writes beside its results, and the only route by
which an `R_T(A_H)` curve enters a run. The manifest declares the tabulated quantity, so that
`R_T` and `r_ν`, which occupy the same column of two files written side by side, cannot be
confused; each curve is resolved through its `temperature_ratio_file` only. `columns` is recorded
as stated by the manifest and is not used to locate a column, since readers take columns by
position. Built by [`read_temperature_ratio_manifest`](@ref).
"""
struct TemperatureRatioManifest
    system::ManifestSystem
    ordinate::String
    abscissa::Vector{String}
    columns::Vector{String}
    curves::Vector{ManifestCurve}
    source::String

    # Concrete argument types, so no `Any`-accepting fallback is generated.
    function TemperatureRatioManifest(
        system::ManifestSystem,
        ordinate::AbstractString,
        abscissa::Vector{String},
        columns::Vector{String},
        curves::Vector{ManifestCurve},
        source::AbstractString,
    )
        return new(system, String(ordinate), abscissa, columns, curves, String(source))
    end
end

function Base.show(io::IO, manifest::TemperatureRatioManifest)
    print(
        io,
        "TemperatureRatioManifest(",
        manifest.system.label,
        ", ",
        length(manifest.curves),
        " curves, ",
        manifest.ordinate,
        " vs ",
        join(manifest.abscissa, ", "),
        ")",
    )
    return nothing
end

# Read one manifest key with its type; a malformed manifest fails here, naming the key.
# `label` is the fully qualified key name used in the message.
function _manifest_value(
    table::AbstractDict,
    key::AbstractString,
    ::Type{T},
    label::AbstractString,
    source::AbstractString,
) where {T}
    haskey(table, key) || throw(ArgumentError("$source has no $label"))
    raw = table[key]
    raw isa T || throw(
        ArgumentError(
            "$source: $label must be $(T), got $(repr(raw)) of type $(typeof(raw))",
        ),
    )
    return raw
end

function _manifest_number(
    table::AbstractDict,
    key::AbstractString,
    label::AbstractString,
    source::AbstractString,
)
    value = _manifest_value(table, key, Real, label, source)
    isfinite(value) || throw(ArgumentError("$source: $label must be finite, got $value"))
    return Float64(value)
end

function _manifest_integer(
    table::AbstractDict,
    key::AbstractString,
    label::AbstractString,
    source::AbstractString,
)
    return Int(_manifest_value(table, key, Int, label, source))
end

"""
    read_temperature_ratio_manifest(path) -> TemperatureRatioManifest

Read and validate the run record of a temperature-ratio extraction. Each check throws an
`ArgumentError` naming the offending key:

- `[run] ordinate` must be `"R_T"` and `[run] abscissa` must be `["A_H"]`.
- Every `[[segmented_curve]]` must carry `label`, `kind`, `temperature_ratio_file` and
  `multiplicity_ratio_pivots_file`, with `kind` drawn from [`MANIFEST_CURVE_KINDS`](@ref).
- Labels must be distinct, since a label selects a curve.
- The two files of a curve must be distinct paths, since they hold different quantities.

`[run] columns` is recorded but not used to locate a column. File paths inside the manifest are
resolved relative to the manifest's directory; absolute paths are taken as written.
"""
function read_temperature_ratio_manifest(path::AbstractString)
    isfile(path) || throw(ArgumentError("no manifest at $path"))
    source = abspath(path)
    document = TOML.parsefile(path)

    for name in ("system", "run")
        haskey(document, name) || throw(ArgumentError("$source has no [$name] table"))
        document[name] isa AbstractDict || throw(
            ArgumentError(
                "$source: [$name] must be a table, got $(typeof(document[name]))",
            ),
        )
    end
    system_table = document["system"]
    run_table = document["run"]

    ordinate = _manifest_value(run_table, "ordinate", String, "[run] ordinate", source)
    ordinate == MANIFEST_ORDINATE || throw(
        ArgumentError("$source: [run] ordinate must be $(repr(MANIFEST_ORDINATE)), got \
 $(repr(ordinate)); this manifest indexes a different quantity"),
    )
    abscissa_raw =
        _manifest_value(run_table, "abscissa", AbstractVector, "[run] abscissa", source)
    abscissa = String[]
    for entry in abscissa_raw
        entry isa AbstractString || throw(
            ArgumentError(
                "$source: [run] abscissa must be a list of strings, got $(repr(entry))",
            ),
        )
        push!(abscissa, String(entry))
    end
    abscissa == MANIFEST_ABSCISSA || throw(
        ArgumentError(
            "$source: [run] abscissa must be $(MANIFEST_ABSCISSA), got $(abscissa)",
        ),
    )

    # Recorded only, never a lookup key; empty when the manifest omits it.
    columns = String[]
    if haskey(run_table, "columns")
        for entry in
            _manifest_value(run_table, "columns", AbstractVector, "[run] columns", source)
            push!(columns, string(entry))
        end
    end

    system = ManifestSystem(
        _manifest_value(system_table, "label", String, "[system] label", source),
        _manifest_value(system_table, "notation", String, "[system] notation", source),
        _manifest_integer(system_table, "target_A", "[system] target_A", source),
        _manifest_integer(system_table, "target_Z", "[system] target_Z", source),
        _manifest_value(system_table, "channel", String, "[system] channel", source),
        _manifest_value(system_table, "reaction", String, "[system] reaction", source),
        _manifest_number(
            system_table,
            "incident_energy_MeV",
            "[system] incident_energy_MeV",
            source,
        ),
        _manifest_integer(system_table, "compound_A", "[system] compound_A", source),
        _manifest_integer(system_table, "compound_Z", "[system] compound_Z", source),
    )

    haskey(document, "segmented_curve") || throw(
        ArgumentError("$source declares no [[segmented_curve]]; there is no curve to take"),
    )
    entries = document["segmented_curve"]
    entries isa AbstractVector || throw(
        ArgumentError(
            "$source: [[segmented_curve]] must be an array of tables, got $(typeof(entries))",
        ),
    )
    isempty(entries) && throw(
        ArgumentError("$source declares no [[segmented_curve]]; there is no curve to take"),
    )

    curves = ManifestCurve[]
    for (index, entry) in pairs(entries)
        entry isa AbstractDict ||
            throw(ArgumentError("$source: [[segmented_curve]] $index must be a table"))
        where_ = "[[segmented_curve]] $index"
        label = _manifest_value(entry, "label", String, "$where_ label", source)
        isempty(label) && throw(ArgumentError("$source: $where_ label must not be empty"))
        kind = _manifest_value(entry, "kind", String, "$where_ kind", source)
        kind in MANIFEST_CURVE_KINDS ||
            throw(ArgumentError("$source: $where_ kind must be one of \
 $(join(map(repr, MANIFEST_CURVE_KINDS), ", ")), got $(repr(kind))"))
        temperature_ratio_file = _manifest_value(
            entry,
            "temperature_ratio_file",
            String,
            "$where_ temperature_ratio_file",
            source,
        )
        pivots_file = _manifest_value(
            entry,
            "multiplicity_ratio_pivots_file",
            String,
            "$where_ multiplicity_ratio_pivots_file",
            source,
        )
        temperature_ratio_file == pivots_file && throw(
            ArgumentError(
                "$source: $where_ names one file as both temperature_ratio_file and \
 multiplicity_ratio_pivots_file, $(repr(pivots_file)); those hold R_T and r_nu \
 and cannot be the same table",
            ),
        )
        push!(curves, ManifestCurve(label, kind, temperature_ratio_file, pivots_file))
    end

    labels = [curve.label for curve in curves]
    allunique(labels) || throw(
        ArgumentError(
            "$source repeats a [[segmented_curve]] label, and a label is what selects a \
 curve: $(join(map(repr, labels), ", "))",
        ),
    )

    return TemperatureRatioManifest(system, ordinate, abscissa, columns, curves, source)
end

"""
    staged_manifest(directory) -> String

The path of the single `manifest_<run>.toml` in `directory`, the staged copy of a
temperature-ratio run's results directory. The run token stays in the file name, so the staged
run is traceable without the configuration naming it. A directory holding no manifest, or more
than one, is refused with an error naming the directory.
"""
function staged_manifest(directory::AbstractString)
    isdir(directory) ||
        throw(ArgumentError("the staged temperature-ratio run does not exist: $directory"))
    names = sort(
        filter(
            name -> startswith(name, "manifest_") && endswith(name, ".toml"),
            readdir(directory),
        ),
    )
    isempty(names) && throw(
        ArgumentError(
            "$directory holds no manifest_<run>.toml; stage the results directory of a temperature-ratio run here, manifest and tables together",
        ),
    )
    length(names) == 1 || throw(
        ArgumentError(
            "$directory holds $(length(names)) manifests and nothing says which run to take: $(join(names, ", ")); keep one staged run per system",
        ),
    )
    return joinpath(directory, only(names))
end

"""
    curve_labels(manifest) -> Vector{String}

The label of every `[[segmented_curve]]` the manifest lists, in the order it lists them.
"""
curve_labels(manifest::TemperatureRatioManifest) =
    [curve.label for curve in manifest.curves]

"""
    manifest_curve(manifest, label) -> ManifestCurve

The curve the manifest lists under `label`, or an error naming the labels it does list.
"""
function manifest_curve(manifest::TemperatureRatioManifest, label::AbstractString)
    for curve in manifest.curves
        curve.label == label && return curve
    end
    throw(
        ArgumentError(
            "$(manifest.source) lists no [[segmented_curve]] labelled $(repr(label)); \
 it lists $(join(map(repr, curve_labels(manifest)), ", "))",
        ),
    )
end

"""
    temperature_ratio_path(manifest, label) -> String

The file holding `R_T(A_H)` for the curve labelled `label`, resolved against the directory the
manifest sits in. This is the only route from a manifest to a file. No counterpart exists for
`multiplicity_ratio_pivots_file`, which holds the multiplicity ratio `r_ν` in the same column
position.
"""
function temperature_ratio_path(manifest::TemperatureRatioManifest, label::AbstractString)
    curve = manifest_curve(manifest, label)
    relative = curve.temperature_ratio_file
    path = isabspath(relative) ? relative : joinpath(dirname(manifest.source), relative)
    isfile(path) || throw(
        ArgumentError("$(manifest.source) lists $(repr(relative)) for the curve labelled \
 $(repr(label)), and no file is there: $path"),
    )
    return path
end
