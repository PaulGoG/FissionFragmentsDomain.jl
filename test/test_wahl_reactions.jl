# The charge distribution a system gets when no evaluated table is supplied, resolved in layers:
# this reaction's own fitted parameters, else the systematics, else the conventional means. Each
# layer has to say which it is, because which one answered changes how much the result is worth.

@testset "per-reaction parameters" begin
    # Table A covers exactly the four reactions the 1988 evaluation was made for.
    @test length(WAHL_PER_REACTION) == 4
    @test Set(p.label for p in values(WAHL_PER_REACTION)) ==
          Set(["U235T", "U233T", "PU239T", "CF252S"])

    for ((Z_F, A_F), parameters) in WAHL_PER_REACTION
        @test parameters.Z_F == Z_F
        @test parameters.A_F == A_F
        # Every parameter that enters a Gaussian has to be usable as one.
        @test parameters.σZ > 0
        @test parameters.σZ50 > 0
        # The shell region is narrower than the peak, which is the Z = 50 closure showing.
        @test parameters.σZ50 < parameters.σZ
        # The even-odd factors are enhancements, so at least unity.
        @test parameters.F_Z >= 1
        @test parameters.F_N >= 1
        # ΔZ is negative for the heavy fragment at A' = 140 and its slope carries it lower.
        @test parameters.ΔZ140 < 0
        @test parameters.ΔZSL < 0
        @test all(in((:σZ50, :F_N, :ΔA_Z, :ΔZmax)), parameters.assumed)
        # Every entry says which entrance channels its fit covers.
        @test !isempty(parameters.channels)
        @test all(in(CHANNELS), parameters.channels)
    end

    # U235T is the reaction with the most data, so nothing in it was assumed.
    @test isempty(WAHL_PER_REACTION[(92, 236)].assumed)
    # CF252S is the sparsest, and the table marks four of its parameters as assumed.
    @test length(WAHL_PER_REACTION[(98, 252)].assumed) == 4

    # Matched on the fissioning nucleus and the entrance channel. The resonance channel differs
    # from thermal by the incident energy alone, so it takes the thermal fit.
    thermal = neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth")
    resonance = neutron_induced_fission(Nuclide(92, 235), 5.8013e-4, "nres")
    @test wahl_reaction_parameters(thermal).label == "U235T"
    @test wahl_reaction_parameters(resonance).label == "U235T"
    @test wahl_reaction_parameters(spontaneous_fission(Nuclide(98, 252))).label == "CF252S"

    # The nucleus alone is not enough. 240-Pu(sf) reaches 239-Pu(nth,f)'s compound nucleus
    # 6.5 MeV lower, where the fitted parameters do not hold, so it gets none and falls through.
    @test wahl_reaction_parameters(spontaneous_fission(Nuclide(94, 240))) === nothing
    # Nor does a fast-neutron channel on a nucleus the table covers.
    @test wahl_reaction_parameters(
        neutron_induced_fission(Nuclide(92, 235), 2.0, "nfast"),
    ) === nothing
    # And 252-Cf's fit is spontaneous, so a neutron-induced route to the same nucleus is not it.
    @test wahl_reaction_parameters(
        neutron_induced_fission(Nuclide(98, 251), 2.53e-8, "nth"),
    ) === nothing

    # A reaction the evaluation does not cover has none.
    @test wahl_reaction_parameters(spontaneous_fission(Nuclide(100, 256))) === nothing
end

@testset "the charge distribution resolves in layers" begin
    table = mass_table()
    heavy(system) = begin
        lower = ceil(Int, system.compound.A / 2)
        lower:(lower + 25)
    end

    # A covered reaction takes its own parameters, and says so with its citation.
    covered = neutron_induced_fission(Nuclide(92, 233), 2.53e-8, "nth")
    first_layer = build_charge_distribution(table, covered, heavy(covered))
    @test occursin("per-reaction parameters for U233T", first_layer.source)
    @test occursin("39, 1 (1988)", first_layer.source)

    # A nucleus inside the fitted range of the systematics but with no entry of its own falls
    # through to the second layer.
    uncovered = neutron_induced_fission(Nuclide(90, 232), 1.5, "nfast")
    @test wahl_reaction_parameters(uncovered) === nothing
    @test is_wahl_applicable(table, uncovered)
    second_layer = build_charge_distribution(table, uncovered, heavy(uncovered))
    @test occursin("systematics", second_layer.source)
    @test !occursin("per-reaction", second_layer.source)

    # Outside the fitted range entirely, only the conventional means are left.
    outside = spontaneous_fission(Nuclide(100, 256))
    @test !is_wahl_applicable(table, outside)
    third_layer = build_charge_distribution(table, outside, heavy(outside))
    @test third_layer.source == "mean values"
    @test charge_polarization(third_layer, 130) == DEFAULT_CHARGE_POLARIZATION
    @test charge_dispersion(third_layer, 130) == DEFAULT_CHARGE_DISPERSION

    # Every layer returns the same type, so nothing downstream has to know which answered.
    for d in (first_layer, second_layer, third_layer)
        @test d isa ChargeDistribution
        @test charge_dispersion(d, 130) > 0
    end
end

guarded(
    CHARGE_DISTRIBUTION_AVAILABLE,
    "the first layer beats the second on every evaluated table",
) do
    @testset "the first layer beats the second on every evaluated table" begin
        # The point of the layering: where a reaction has its own parameters they are closer to
        # what an evaluation of that reaction gives than a trend across reactions is.
        table = mass_table()
        cases = (
            (
                U235_CHARGE_DISTRIBUTION_FILE,
                neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth"),
            ),
            (CF252_CHARGE_DISTRIBUTION_FILE, spontaneous_fission(Nuclide(98, 252))),
        )
        for (path, system) in cases
            isfile(path) || continue
            evaluated = read_charge_distribution(path)
            masses = sort(collect(keys(evaluated.σ_Z)))
            systematics = wahl_charge_distribution(table, system, masses)
            per_reaction = build_charge_distribution(table, system, masses)
            @test occursin("per-reaction", per_reaction.source)

            mad(d, accessor) =
                sum(abs(accessor(d, A) - accessor(evaluated, A)) for A in masses) /
                length(masses)
            @test mad(per_reaction, charge_polarization) <
                  mad(systematics, charge_polarization)
            @test mad(per_reaction, charge_dispersion) < mad(systematics, charge_dispersion)
        end
    end
end

@testset "Table A places its own regions" begin
    # Footnotes a and c give the ranges the narrow dispersion and the unit even-odd factors apply
    # over. They are part of the parameter set and travel with it, because the widths were fitted
    # on the assumption that these are the regions.
    parameters = WAHL_PER_REACTION[(98, 252)]
    systematics = WahlSystematics(parameters, 0.0)
    @test systematics.shell == parameters.shell
    @test systematics.symmetric == parameters.symmetric

    # Inside the shell region the narrow width, outside it the single average. No slope, and no
    # dependence on eq. (10)'s boundaries.
    @test wahl_dispersion(systematics, 126) == parameters.σZ50
    @test wahl_dispersion(systematics, 129) == parameters.σZ50
    @test wahl_dispersion(systematics, 130) == parameters.σZ
    @test wahl_dispersion(systematics, 160) == parameters.σZ
    # Which is what eq. (10) does not do: it holds the narrow width far past where Table A stops.
    @test wahl_dispersion(WahlSystematics(98, 252, 0.0), 135) !=
          wahl_dispersion(WahlSystematics(98, 252, 0.0), 160)
    @test wahl_dispersion(systematics, 135) == wahl_dispersion(systematics, 160)

    # The even-odd factors are unity through the near-symmetry range and Table A's values outside.
    @test wahl_even_odd_factors(systematics, 126) == (1.0, 1.0)
    @test wahl_even_odd_factors(systematics, 129) == (1.0, 1.0)
    @test wahl_even_odd_factors(systematics, 130) == (parameters.F_Z, parameters.F_N)

    # A systematics evaluation names no regions, so eq. (10) places them as it always has.
    plain = WahlSystematics(98, 252, 0.0)
    @test isempty(plain.shell)
    @test isempty(plain.symmetric)
end
