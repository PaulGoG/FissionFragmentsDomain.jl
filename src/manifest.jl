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
    ManifestDomain

The fragmentation domain an extraction of `R_T(A_H)` was computed on, as its `[domain]` table
records it: the level density model (`level_density_model`, one of
[`LEVEL_DENSITY_MODELS`](@ref)), how the level density ratio entered the inversion
(`ratio_averaging`, one of [`RATIO_AVERAGINGS`](@ref)), the charges retained per mass
(`charges_per_mass`, odd and positive), the charge model (`charge_model`, as
[`charge_model_label`](@ref) spells it), the mass table (`mass_table`, the file name of its
source), the version of this package (`package_version`), whether the Gilbert–Cameron
formula took its deformed branch (`deformed_branch`, required with `"GC"` and `false` otherwise;
see [`GilbertCameron`](@ref)), and whether a charge-resolved
inversion weighted each fragmentation by its mean total excitation (`excitation_weighted`, see
[`ChargeResolved`](@ref); `false` when the table omits it). A consumer that partitions the
excitation energy on a different domain does not invert the same relation.
"""
struct ManifestDomain
    level_density_model::String
    ratio_averaging::String
    charges_per_mass::Int
    charge_model::String
    mass_table::String
    package_version::String
    excitation_weighted::Bool
    deformed_branch::Bool

    function ManifestDomain(
        level_density_model::AbstractString,
        ratio_averaging::AbstractString,
        charges_per_mass::Integer,
        charge_model::AbstractString,
        mass_table::AbstractString,
        package_version::AbstractString,
        excitation_weighted::Bool = false,
        deformed_branch::Bool = false,
    )
        level_density_model in LEVEL_DENSITY_MODELS ||
            throw(ArgumentError("[domain] level_density_model must be one of \
                 $(join(map(repr, LEVEL_DENSITY_MODELS), ", ")), got \
                 $(repr(level_density_model))"))
        ratio_averaging in RATIO_AVERAGINGS || throw(
            ArgumentError(
                "[domain] ratio_averaging must be one of \
                 $(join(map(repr, RATIO_AVERAGINGS), ", ")), got $(repr(ratio_averaging))",
            ),
        )
        (charges_per_mass > 0 && isodd(charges_per_mass)) || throw(
            ArgumentError(
                "[domain] charges_per_mass must be odd and positive, got $charges_per_mass",
            ),
        )
        isempty(charge_model) &&
            throw(ArgumentError("[domain] charge_model must not be empty"))
        isempty(mass_table) && throw(ArgumentError("[domain] mass_table must not be empty"))
        isempty(package_version) &&
            throw(ArgumentError("[domain] package_version must not be empty"))
        deformed_branch &&
            level_density_model != "GC" &&
            throw(
                ArgumentError(
                    "[domain] deformed_branch applies to level_density_model = \"GC\" only, \
                     got $(repr(level_density_model))",
                ),
            )
        excitation_weighted &&
            ratio_averaging != "charge_resolved" &&
            throw(
                ArgumentError(
                    "[domain] excitation_weighted applies to ratio_averaging = \"charge_resolved\" \
                 only, got $(repr(ratio_averaging))",
                ),
            )
        return new(
            String(level_density_model),
            String(ratio_averaging),
            Int(charges_per_mass),
            String(charge_model),
            String(mass_table),
            String(package_version),
            excitation_weighted,
            deformed_branch,
        )
    end
end

"""
    ManifestDomain(model, averaging, domain, charge, masses) -> ManifestDomain

The `[domain]` record of an extraction computed with the level density model `model`, the ratio
averaging `averaging`, on the fragmentation domain `domain` built from the charge model `charge`
and the mass table `masses`, by the version of this package that is loaded.
"""
function ManifestDomain(
    model::LevelDensityModel,
    averaging::RatioAveraging,
    domain::FragmentationDomain,
    charge::ChargeModel,
    masses::MassExcessTable,
)
    return ManifestDomain(
        level_density_label(model),
        ratio_averaging_label(averaging),
        domain.charges_per_mass,
        charge_model_label(charge),
        basename(masses.source),
        string(PACKAGE_VERSION),
        averaging isa ChargeResolved && averaging.excitation !== nothing,
        model isa GilbertCameron && model.deformed_branch,
    )
end

Base.:(==)(a::ManifestDomain, b::ManifestDomain) =
    all(getfield(a, f) == getfield(b, f) for f in fieldnames(ManifestDomain))

Base.hash(domain::ManifestDomain, h::UInt) = hash(
    ntuple(i -> getfield(domain, i), fieldcount(ManifestDomain)),
    hash(ManifestDomain, h),
)

function Base.show(io::IO, domain::ManifestDomain)
    print(
        io,
        "ManifestDomain(",
        domain.level_density_model,
        ", ",
        domain.ratio_averaging,
        ", ",
        domain.charges_per_mass,
        " Z per A, ",
        domain.charge_model,
        ", ",
        domain.mass_table,
        domain.deformed_branch ? ", deformed branch" : "",
        domain.excitation_weighted ? ", excitation-weighted" : "",
        ")",
    )
    return nothing
end

"""
    TemperatureRatioManifest

The run record a temperature-ratio extraction writes beside its results, and the only route by
which an `R_T(A_H)` curve enters a run. The manifest declares the tabulated quantity, so that
`R_T` and `r_ν`, which occupy the same column of two files written side by side, cannot be
confused; each curve is resolved through its `temperature_ratio_file` only. `columns` is recorded
as stated by the manifest and is not used to locate a column, since readers take columns by
position. `domain` holds the optional `[domain]` table, a [`ManifestDomain`](@ref), and is
`nothing` for manifests written before that table existed. Built by
[`read_temperature_ratio_manifest`](@ref).
"""
struct TemperatureRatioManifest
    system::ManifestSystem
    ordinate::String
    abscissa::Vector{String}
    columns::Vector{String}
    curves::Vector{ManifestCurve}
    source::String
    domain::Union{ManifestDomain, Nothing}

    # Concrete argument types, so no `Any`-accepting fallback is generated.
    function TemperatureRatioManifest(
        system::ManifestSystem,
        ordinate::AbstractString,
        abscissa::Vector{String},
        columns::Vector{String},
        curves::Vector{ManifestCurve},
        source::AbstractString,
        domain::Union{ManifestDomain, Nothing},
    )
        return new(
            system,
            String(ordinate),
            abscissa,
            columns,
            curves,
            String(source),
            domain,
        )
    end

    function TemperatureRatioManifest(
        system::ManifestSystem,
        ordinate::AbstractString,
        abscissa::Vector{String},
        columns::Vector{String},
        curves::Vector{ManifestCurve},
        source::AbstractString,
    )
        return TemperatureRatioManifest(
            system,
            ordinate,
            abscissa,
            columns,
            curves,
            source,
            nothing,
        )
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

# Required with Gilbert–Cameron, whose record would otherwise read as one branch or the other.
function _manifest_deformed_branch(table::AbstractDict, source::AbstractString)
    haskey(table, "deformed_branch") && return _manifest_value(
        table,
        "deformed_branch",
        Bool,
        "[domain] deformed_branch",
        source,
    )
    get(table, "level_density_model", nothing) == "GC" && throw(
        ArgumentError(
            "$source: [domain] deformed_branch is required with level_density_model = \"GC\"",
        ),
    )
    return false
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
- `[domain]` is optional; when present it must be a table carrying `level_density_model`, one
  of [`LEVEL_DENSITY_MODELS`](@ref), `ratio_averaging`, one of [`RATIO_AVERAGINGS`](@ref), an
  odd positive integer `charges_per_mass`, and non-empty `charge_model`, `mass_table` and
  `package_version`; `excitation_weighted`, a Boolean, is optional and admitted only with
  `ratio_averaging = "charge_resolved"`; `deformed_branch`, a Boolean, is required with
  `level_density_model = "GC"` and must be `false` or absent otherwise. Without the table the manifest's `domain` is `nothing`.

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

    # Optional: manifests written before the table existed carry no domain record.
    domain = nothing
    if haskey(document, "domain")
        domain_table = document["domain"]
        domain_table isa AbstractDict || throw(
            ArgumentError("$source: [domain] must be a table, got $(typeof(domain_table))"),
        )
        fields = (
            _manifest_value(
                domain_table,
                "level_density_model",
                String,
                "[domain] level_density_model",
                source,
            ),
            _manifest_value(
                domain_table,
                "ratio_averaging",
                String,
                "[domain] ratio_averaging",
                source,
            ),
            _manifest_integer(
                domain_table,
                "charges_per_mass",
                "[domain] charges_per_mass",
                source,
            ),
            _manifest_value(
                domain_table,
                "charge_model",
                String,
                "[domain] charge_model",
                source,
            ),
            _manifest_value(
                domain_table,
                "mass_table",
                String,
                "[domain] mass_table",
                source,
            ),
            _manifest_value(
                domain_table,
                "package_version",
                String,
                "[domain] package_version",
                source,
            ),
            haskey(domain_table, "excitation_weighted") ?
            _manifest_value(
                domain_table,
                "excitation_weighted",
                Bool,
                "[domain] excitation_weighted",
                source,
            ) : false,
            _manifest_deformed_branch(domain_table, source),
        )
        domain = try
            ManifestDomain(fields...)
        catch err
            err isa ArgumentError || rethrow()
            throw(ArgumentError("$source: " * err.msg))
        end
    end

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

    return TemperatureRatioManifest(
        system,
        ordinate,
        abscissa,
        columns,
        curves,
        source,
        domain,
    )
end

# The run record's tables in reading order; keys within a table alphabetically.
const _MANIFEST_TABLE_ORDER = ("system", "run", "domain", "segmented_curve")
_manifest_key_order(key::AbstractString) =
    (something(findfirst(==(key), _MANIFEST_TABLE_ORDER), 0), key)

"""
    write_temperature_ratio_manifest(path, manifest::TemperatureRatioManifest) -> String

Write `manifest` as the TOML run record that [`read_temperature_ratio_manifest`](@ref)
reads: `[system]`, `[run]`, `[domain]` when present, then one `[[segmented_curve]]` per
curve, file paths as the curves hold them. Refuses to overwrite an existing file, so a run
record is never lost; returns `path`. What it writes reads back equal, `source` aside.
"""
function write_temperature_ratio_manifest(
    path::AbstractString,
    manifest::TemperatureRatioManifest,
)
    isfile(path) && throw(ArgumentError("$path exists; a run record is not overwritten"))
    system = manifest.system
    run_table =
        Dict{String, Any}("ordinate" => manifest.ordinate, "abscissa" => manifest.abscissa)
    isempty(manifest.columns) || (run_table["columns"] = manifest.columns)
    document = Dict{String, Any}(
        "system" => Dict{String, Any}(
            "label" => system.label,
            "notation" => system.notation,
            "target_A" => system.target_A,
            "target_Z" => system.target_Z,
            "channel" => system.channel,
            "reaction" => system.reaction,
            "incident_energy_MeV" => system.incident_energy_MeV,
            "compound_A" => system.compound_A,
            "compound_Z" => system.compound_Z,
        ),
        "run" => run_table,
        "segmented_curve" => [
            Dict{String, Any}(
                "label" => curve.label,
                "kind" => curve.kind,
                "temperature_ratio_file" => curve.temperature_ratio_file,
                "multiplicity_ratio_pivots_file" =>
                    curve.multiplicity_ratio_pivots_file,
            ) for curve in manifest.curves
        ],
    )
    domain = manifest.domain
    if domain !== nothing
        document["domain"] = Dict{String, Any}(
            "level_density_model" => domain.level_density_model,
            "ratio_averaging" => domain.ratio_averaging,
            "charges_per_mass" => domain.charges_per_mass,
            "charge_model" => domain.charge_model,
            "mass_table" => domain.mass_table,
            "package_version" => domain.package_version,
            "excitation_weighted" => domain.excitation_weighted,
            "deformed_branch" => domain.deformed_branch,
        )
    end
    open(io -> TOML.print(io, document; sorted = true, by = _manifest_key_order), path, "w")
    return path
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
