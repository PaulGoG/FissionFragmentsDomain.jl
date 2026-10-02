# The manifest is what makes the quantity a matter of record rather than of staging. These
# exercise every refusal, because a manifest that is read leniently is worth no more than the
# hand-staged path it replaced.

"Write a manifest with `replacements` applied to the fixture, and the two tables it names."
function _manifest(directory::AbstractString, replacements::Pair{String, String}...)
    text = """
    [system]
    label = "U235_nth"
    notation = "²³⁵U(nth,f)"
    target_A = 235
    target_Z = 92
    channel = "nth"
    reaction = "n,f"
    incident_energy_MeV = 2.53e-8
    compound_A = 236
    compound_Z = 92

    [run]
    ordinate = "R_T"
    abscissa = ["A_H"]
    columns = ["A_H", "R_T", "R_T_uncertainty"]

    [[segmented_curve]]
    label = "Goeoek_2018"
    kind = "dataset"
    temperature_ratio_file = "R_T_vs_A_H_segmented_Goeoek_2018_run.csv"
    multiplicity_ratio_pivots_file = "r_nu_vs_A_H_pivots_Goeoek_2018_run.csv"

    [[segmented_curve]]
    label = "systematic_trend"
    kind = "systematic_trend"
    temperature_ratio_file = "R_T_vs_A_H_segmented_systematic_trend_run.csv"
    multiplicity_ratio_pivots_file = "r_nu_vs_A_H_pivots_systematic_trend_run.csv"
    """
    for (from, to) in replacements
        occursin(from, text) || error("the fixture has no $(repr(from)) to replace")
        text = replace(text, from => to)
    end
    # R_T rises through unity; r_nu sits at one half. Both are plausible floats in the same
    # column, which is the whole reason the manifest exists.
    for label in ("Goeoek_2018", "systematic_trend")
        write(
            joinpath(directory, "R_T_vs_A_H_segmented_$(label)_run.csv"),
            "A_H,R_T,R_T_uncertainty\n118,1.0,0.0\n140,1.2,0.01\n160,1.4,0.02\n",
        )
        write(
            joinpath(directory, "r_nu_vs_A_H_pivots_$(label)_run.csv"),
            "A_H,r_nu,r_nu_uncertainty\n118,0.5,0.0\n140,0.5,0.01\n160,0.5,0.02\n",
        )
    end
    path = joinpath(directory, "manifest_run.toml")
    write(path, text)
    return path
end

function _manifest_rejects(
    directory::AbstractString,
    key::AbstractString,
    replacements::Pair{String, String}...,
)
    thrown = try
        read_temperature_ratio_manifest(_manifest(directory, replacements...))
        nothing
    catch err
        err
    end
    thrown isa ArgumentError || return false
    return occursin(key, thrown.msg)
end

@testset "reading a temperature ratio manifest" begin
    mktempdir() do directory
        manifest = read_temperature_ratio_manifest(_manifest(directory))
        @test manifest isa TemperatureRatioManifest
        @test manifest.ordinate == MANIFEST_ORDINATE
        @test manifest.abscissa == MANIFEST_ABSCISSA
        @test curve_labels(manifest) == ["Goeoek_2018", "systematic_trend"]
        @test manifest.system.label == "U235_nth"
        @test manifest.system.compound_A == 236
        @test manifest.system.incident_energy_MeV ≈ 2.53e-8 rtol = RTOL

        trend = manifest_curve(manifest, "systematic_trend")
        @test trend.kind == "systematic_trend"
        @test manifest_curve(manifest, "Goeoek_2018").kind == "dataset"

        # The resolved path is the temperature-ratio table and never the pivots table, and it is
        # taken relative to the directory the manifest sits in.
        resolved = temperature_ratio_path(manifest, "systematic_trend")
        @test basename(resolved) == trend.temperature_ratio_file
        @test basename(resolved) != trend.multiplicity_ratio_pivots_file
        @test dirname(resolved) == directory
        @test isfile(resolved)

        # And what it resolves to really is R_T rather than r_nu: the curve rises through unity,
        # where the pivots table beside it sits flat at one half.
        curve = read_segmented_curve(resolved)
        @test curve(118) ≈ 1.0 rtol = RTOL
        @test curve(160) ≈ 1.4 rtol = RTOL
        @test curve.source == resolved

        @test_throws ArgumentError manifest_curve(manifest, "absent")
        @test_throws ArgumentError temperature_ratio_path(manifest, "absent")
        # The message names the labels on offer rather than only the one that was asked for.
        message = try
            manifest_curve(manifest, "absent")
        catch err
            err.msg
        end
        @test occursin("systematic_trend", message)

        @test_throws ArgumentError read_temperature_ratio_manifest(
            joinpath(directory, "absent.toml"),
        )
    end
end

@testset "a manifest this package cannot honour is refused" begin
    mktempdir() do directory
        # A manifest indexing another quantity belongs to another consumer. This is the check
        # that makes staging the multiplicity ratio impossible rather than merely inadvisable.
        @test _manifest_rejects(
            directory,
            "ordinate",
            "ordinate = \"R_T\"" => "ordinate = \"r_nu\"",
        )
        @test _manifest_rejects(
            directory,
            "abscissa",
            "abscissa = [\"A_H\"]" => "abscissa = [\"A\"]",
        )
        @test _manifest_rejects(
            directory,
            "kind",
            "kind = \"systematic_trend\"" => "kind = \"guess\"",
        )
        # A label is what selects a curve, so two of them cannot share one.
        @test _manifest_rejects(
            directory,
            "label",
            "label = \"Goeoek_2018\"" => "label = \"systematic_trend\"",
        )
        # One path under both keys means the producer lost the distinction before this package
        # saw it, and there is then no way to tell which quantity the file holds.
        @test _manifest_rejects(
            directory,
            "multiplicity_ratio_pivots_file",
            "multiplicity_ratio_pivots_file = \"r_nu_vs_A_H_pivots_systematic_trend_run.csv\"" => "multiplicity_ratio_pivots_file = \"R_T_vs_A_H_segmented_systematic_trend_run.csv\"",
        )
        @test _manifest_rejects(directory, "[run] ordinate", "ordinate = \"R_T\"\n" => "")
        @test _manifest_rejects(directory, "[system] label", "label = \"U235_nth\"\n" => "")
        @test _manifest_rejects(
            directory,
            "[system] target_A",
            "target_A = 235" => "target_A = \"235\"",
        )
        @test _manifest_rejects(
            directory,
            "segmented_curve",
            "[[segmented_curve]]\nlabel = \"Goeoek_2018\"\nkind = \"dataset\"\ntemperature_ratio_file = \"R_T_vs_A_H_segmented_Goeoek_2018_run.csv\"\nmultiplicity_ratio_pivots_file = \"r_nu_vs_A_H_pivots_Goeoek_2018_run.csv\"\n\n[[segmented_curve]]\nlabel = \"systematic_trend\"\nkind = \"systematic_trend\"\ntemperature_ratio_file = \"R_T_vs_A_H_segmented_systematic_trend_run.csv\"\nmultiplicity_ratio_pivots_file = \"r_nu_vs_A_H_pivots_systematic_trend_run.csv\"\n" => "",
        )

        # A listed file that is not there is caught when the curve is resolved, not silently.
        manifest = read_temperature_ratio_manifest(
            _manifest(
                directory,
                "temperature_ratio_file = \"R_T_vs_A_H_segmented_systematic_trend_run.csv\"" => "temperature_ratio_file = \"absent.csv\"",
            ),
        )
        @test_throws ArgumentError temperature_ratio_path(manifest, "systematic_trend")
    end
end

@testset "finding the staged run" begin
    mktempdir() do directory
        @test_throws ArgumentError staged_manifest(joinpath(directory, "absent"))

        empty_directory = joinpath(directory, "empty")
        mkpath(empty_directory)
        @test_throws ArgumentError staged_manifest(empty_directory)

        one = joinpath(directory, "one")
        mkpath(one)
        _manifest(one)
        @test staged_manifest(one) == joinpath(one, "manifest_run.toml")

        # Two staged runs in one place, and nothing says which to take.
        cp(joinpath(one, "manifest_run.toml"), joinpath(one, "manifest_other.toml"))
        @test_throws ArgumentError staged_manifest(one)
    end
end

@testset "a manifest written from the shared record reads back" begin
    @test SYSTEMATIC_TREND_LABEL in MANIFEST_CURVE_KINDS
    @test !any(isspace, SYSTEMATIC_TREND_LABEL)

    mktempdir() do directory
        system = neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth")
        write(joinpath(directory, "R_T.csv"), "A_H,R_T\n118,1.0\n160,1.4\n")
        write(joinpath(directory, "r_nu.csv"), "A_H,r_nu\n118,0.5\n160,0.4\n")
        document = Dict{String, Any}(
            "system" => system_record(system),
            "run" => Dict{String, Any}(
                "ordinate" => MANIFEST_ORDINATE,
                "abscissa" => MANIFEST_ABSCISSA,
                "columns" => ["A_H", "R_T", "R_T_uncertainty"],
            ),
            "segmented_curve" => [
                Dict{String, Any}(
                    "label" => SYSTEMATIC_TREND_LABEL,
                    "kind" => "systematic_trend",
                    "temperature_ratio_file" => "R_T.csv",
                    "multiplicity_ratio_pivots_file" => "r_nu.csv",
                ),
            ],
        )
        path = joinpath(directory, "manifest_run.toml")
        open(io -> TOML.print(io, document), path, "w")
        manifest = read_temperature_ratio_manifest(path)
        @test manifest.system.label == system_label(system)
        @test curve_labels(manifest) == [SYSTEMATIC_TREND_LABEL]
        @test read_segmented_curve(temperature_ratio_path(manifest, SYSTEMATIC_TREND_LABEL))(
            160,
        ) == 1.4
    end
end

@testset "the archive accession of a dataset curve is recorded" begin
    @test ManifestCurve("x", "dataset", "a.csv", "b.csv").accession == ""
    for accession in ("14369003", "V0101002", "309690021")
        @test ManifestCurve("x", "dataset", "a.csv", "b.csv"; accession).accession ==
              accession
    end
    for accession in ("1436900", "1436900345", "14369 003", "v0101002")
        @test_throws ArgumentError ManifestCurve(
            "x",
            "dataset",
            "a.csv",
            "b.csv";
            accession,
        )
    end
    @test_throws ArgumentError ManifestCurve(
        SYSTEMATIC_TREND_LABEL,
        "systematic_trend",
        "a.csv",
        "b.csv";
        accession = "14369003",
    )

    mktempdir() do directory
        # Absent, the key reads as empty and is not written back.
        base = read_temperature_ratio_manifest(_manifest(directory))
        @test all(isempty(curve.accession) for curve in base.curves)
        plain = write_temperature_ratio_manifest(joinpath(directory, "plain.toml"), base)
        @test !occursin("accession", read(plain, String))
    end

    mktempdir() do directory
        path = _manifest(
            directory,
            "kind = \"dataset\"\n" => "kind = \"dataset\"\naccession = \"23444004\"\n",
        )
        manifest = read_temperature_ratio_manifest(path)
        @test manifest_curve(manifest, "Goeoek_2018").accession == "23444004"
        @test manifest_curve(manifest, SYSTEMATIC_TREND_LABEL).accession == ""
        written = write_temperature_ratio_manifest(
            joinpath(directory, "manifest_written.toml"),
            manifest,
        )
        back = read_temperature_ratio_manifest(written)
        @test [curve.accession for curve in back.curves] == ["23444004", ""]
        @test curve_labels(back) == curve_labels(manifest)
    end

    mktempdir() do directory
        @test _manifest_rejects(
            directory,
            "accession",
            "kind = \"dataset\"\n" => "kind = \"dataset\"\naccession = \"2344\"\n",
        )
    end
    mktempdir() do directory
        @test _manifest_rejects(
            directory,
            "accession",
            "kind = \"systematic_trend\"\n" => "kind = \"systematic_trend\"\naccession = \"23444004\"\n",
        )
    end
    mktempdir() do directory
        @test _manifest_rejects(
            directory,
            "accession",
            "kind = \"dataset\"\n" => "kind = \"dataset\"\naccession = 23444004\n",
        )
    end
end

@testset "the domain record round-trips" begin
    mktempdir() do directory
        system = neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth")
        masses = mass_table()
        model = BackShiftedFermiGas(masses)
        charge = charge_model(masses, system)
        domain = fragmentation_domain(system, charge, 118:160)
        record = ManifestDomain(model, ChargeResolved(), domain, charge, masses)
        @test record.level_density_model == "BSFG"
        @test record.ratio_averaging == "charge_resolved"
        @test record.charges_per_mass == 5
        @test record.mass_table == basename(String(AME2020_MASS_EXCESS_FILE))
        @test occursin("Wahl1988", record.charge_model)

        base = read_temperature_ratio_manifest(_manifest(directory))
        manifest = TemperatureRatioManifest(
            base.system,
            base.ordinate,
            base.abscissa,
            base.columns,
            base.curves,
            base.source,
            record,
        )
        path = write_temperature_ratio_manifest(
            joinpath(directory, "manifest_written.toml"),
            manifest,
        )
        back = read_temperature_ratio_manifest(path)
        for field in fieldnames(ManifestSystem)
            @test getfield(back.system, field) == getfield(base.system, field)
        end
        @test back.ordinate == base.ordinate
        @test back.abscissa == base.abscissa
        @test back.columns == base.columns
        @test [curve.label for curve in back.curves] == [curve.label for curve in base.curves]
        @test [curve.kind for curve in back.curves] == [curve.kind for curve in base.curves]
        @test [curve.temperature_ratio_file for curve in back.curves] == [curve.temperature_ratio_file for curve in base.curves]
        @test [curve.multiplicity_ratio_pivots_file for curve in back.curves] == [curve.multiplicity_ratio_pivots_file for curve in base.curves]
        @test back.domain == record

        # A run record is never overwritten.
        @test_throws ArgumentError write_temperature_ratio_manifest(path, manifest)

        for averaging in (ChargeResolved(), RatioOfMeans(), MeanOfRatios())
            @test ratio_averaging(ratio_averaging_label(averaging)) isa typeof(averaging)
        end
        @test_throws ArgumentError ratio_averaging("median")
    end
end

@testset "a domain record this package cannot honour is refused" begin
    mktempdir() do directory
        block = "[domain]\nlevel_density_model = \"BSFG\"\nratio_averaging = \"charge_resolved\"\ncharges_per_mass = 5\ncharge_model = \"Wahl1988(U235T, Z_F = 92, A_F = 236)\"\nmass_table = \"mass_excess_ame2020.dat\"\npackage_version = \"0.1.3\"\n\n[run]"
        altered(from, to) = "[run]" => replace(block, from => to)

        @test read_temperature_ratio_manifest(_manifest(directory, "[run]" => block)).domain isa
              ManifestDomain
        @test read_temperature_ratio_manifest(_manifest(directory)).domain === nothing

        @test _manifest_rejects(
            directory,
            "ratio_averaging",
            altered(
                "ratio_averaging = \"charge_resolved\"",
                "ratio_averaging = \"median\"",
            ),
        )
        @test _manifest_rejects(
            directory,
            "level_density_model",
            altered("level_density_model = \"BSFG\"", "level_density_model = \"Fermi\""),
        )
        @test _manifest_rejects(
            directory,
            "charges_per_mass",
            altered("charges_per_mass = 5", "charges_per_mass = 4"),
        )
        @test _manifest_rejects(
            directory,
            "charges_per_mass",
            altered("charges_per_mass = 5", "charges_per_mass = \"5\""),
        )
        @test _manifest_rejects(
            directory,
            "mass_table",
            altered("mass_table = \"mass_excess_ame2020.dat\"\n", ""),
        )
    end
end

@testset "the excitation weighting is recorded" begin
    masses = mass_table()
    system = neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth")
    model = BackShiftedFermiGas(masses)
    charge = charge_model(masses, system)
    domain = fragmentation_domain(system, charge, 118:160)
    plain = ManifestDomain(model, ChargeResolved(), domain, charge, masses)
    @test !plain.excitation_weighted
    excitation = mean_total_excitation(masses, domain, Dict(A => 170.0 for A in 118:160))
    weighted = ManifestDomain(model, ChargeResolved(excitation), domain, charge, masses)
    @test weighted.excitation_weighted
    @test weighted != plain
    @test_throws ArgumentError ManifestDomain(
        "BSFG",
        "ratio_of_means",
        5,
        "x",
        "y",
        "0.1.3",
        true,
    )

    mktempdir() do directory
        base = read_temperature_ratio_manifest(_manifest(directory))
        manifest = TemperatureRatioManifest(
            base.system,
            base.ordinate,
            base.abscissa,
            base.columns,
            base.curves,
            base.source,
            weighted,
        )
        path =
            write_temperature_ratio_manifest(joinpath(directory, "weighted.toml"), manifest)
        @test read_temperature_ratio_manifest(path).domain == weighted
        # The record reads in the order a person reads it.
        tables = filter(line -> startswith(line, "["), readlines(path))
        @test first(tables) == "[system]"
        @test findfirst(==("[domain]"), tables) <
              findfirst(==("[[segmented_curve]]"), tables)
    end
end

@testset "the Gilbert–Cameron branch is recorded" begin
    @test_throws ArgumentError ManifestDomain(
        "BSFG",
        "charge_resolved",
        5,
        "x",
        "y",
        "0.1.5",
        false,
        true,
    )
    gc = ManifestDomain("GC", "charge_resolved", 5, "x", "y", "0.1.5", false, true)
    @test gc.deformed_branch
    @test gc != ManifestDomain("GC", "charge_resolved", 5, "x", "y", "0.1.5", false, false)

    mktempdir() do directory
        base = read_temperature_ratio_manifest(_manifest(directory))
        manifest = TemperatureRatioManifest(
            base.system,
            base.ordinate,
            base.abscissa,
            base.columns,
            base.curves,
            base.source,
            gc,
        )
        path = write_temperature_ratio_manifest(joinpath(directory, "gc.toml"), manifest)
        @test read_temperature_ratio_manifest(path).domain == gc

        # A Gilbert–Cameron record that does not say which branch it took is refused.
        write(path, replace(read(path, String), r"deformed_branch = \w+\n" => ""))
        thrown = try
            read_temperature_ratio_manifest(path)
            nothing
        catch err
            err
        end
        @test thrown isa ArgumentError
        @test occursin("[domain] deformed_branch", thrown.msg)
    end
end
