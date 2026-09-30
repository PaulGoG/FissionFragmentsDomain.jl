# Wahl (1988) reproduced from its own tables. The calculated columns of Tables I–IV were produced
# by the model with the parameters of Table A, so a faithful implementation has to return them
# from A and ν̄_A alone, row by row, at the printed precision.

const WAHL_1988_TABLES = Dict(
    label => begin
        table = read_delimited_table(
            joinpath(@__DIR__, "references", "wahl1988", label * ".csv"),
            TableSpec(
                [:A, :NU, :Zp, :AVEZ_calc, :AVEZ_expt, :RMS_calc, :RMS_expt];
                delimiter = ',',
                skip = 1,
            ),
        )
        (
            A = integer_column(table, :A),
            NU = column(table, :NU, Float64),
            Zp = column(table, :Zp, Float64),
            AVEZ = column(table, :AVEZ_calc, Float64),
            RMS = column(table, :RMS_calc, Float64),
            AVEZ_expt = column(table, :AVEZ_expt, Float64),
            RMS_expt = column(table, :RMS_expt, Float64),
        )
    end for label in ("U235T", "U233T", "PU239T", "CF252S")
)

wahl_1988_model(label) = only(m for m in values(WAHL_1988) if m.label == label)

@testset "Table A as printed, and as used" begin
    @test Set(keys(WAHL_1988_TABLE_A)) == Set(["U235T", "U233T", "PU239T", "CF252S"])
    @test WAHL_1988_TABLE_A["U235T"].ΔZ140 == (-0.511, 0.005)
    @test WAHL_1988_TABLE_A["CF252S"].F_N == (1.00, nothing)
    @test WAHL_1988_TABLE_A["PU239T"].σZ50 == (0.35, nothing)

    # The values the tables were computed with lie within the rounding of the printed ones: half
    # a unit of the last printed digit.
    for model in values(WAHL_1988), key in keys(WAHL_1988_TABLE_A[model.label])
        printed = first(getproperty(WAHL_1988_TABLE_A[model.label], key))
        digits =
            something(findfirst(d -> round(printed; digits = d) == printed, 0:4), 5) - 1
        @test abs(getproperty(model, key) - printed) <= 0.5 * 10.0^(-digits) + 1e-12
    end
end

@testset "Wahl (1988), Tables I–IV reproduced" begin
    for (label, table) in WAHL_1988_TABLES
        model = wahl_1988_model(label)
        Zp_error = 0.0
        moment_error = 0.0
        for (i, A) in enumerate(table.A)
            A′ = A + table.NU[i]
            Zp_error = max(Zp_error, abs(most_probable_charge(model, A′) - table.Zp[i]))
            Z̄, RMS = charge_moments(fractional_independent_yields(model, A, table.NU[i])...)
            moment_error =
                max(moment_error, abs(Z̄ - table.AVEZ[i]), abs(RMS - table.RMS[i]))
        end
        # Printed to three decimals: 0.0005 by rounding, and the few masses of U235T at 0.0017
        # are the arithmetic of the original computation, which no parameter moves.
        @test Zp_error < 0.001
        @test moment_error < 0.002
    end

    # With the printed Table A the slopes are rounded, and the extreme masses show it.
    printed = Wahl1988("U233T"; printed = true)
    table = WAHL_1988_TABLES["U233T"]
    departure = maximum(
        abs(most_probable_charge(printed, A + ν) - Zp) for
        (A, ν, Zp) in zip(table.A, table.NU, table.Zp)
    )
    @test departure > 0.01
end

@testset "the near-symmetry regions follow A'" begin
    # Footnote a of Table A lists the product masses that take σ̄₅₀. The tables apply it on the
    # steep branch of ΔZ in A', which for U235T excludes A = 109 although the footnote prints
    # 105–109: its A'_Hc = 124.996 lies below A'_m = 125.117.
    used = Dict(
        "U235T" => [105:108; 125:129],
        "U233T" => [104:108; 123:128],
        "PU239T" => [109:113; 124:129],
        "CF252S" => collect(119:129),
    )
    for (label, masses) in used
        model = wahl_1988_model(label)
        table = WAHL_1988_TABLES[label]
        narrow = [
            A for (A, ν) in zip(table.A, table.NU) if
            charge_dispersion(model, A + ν) == model.σZ50
        ]
        @test narrow == masses
        # F = 1 from the junction down to symmetry, on both sides.
        unity = [
            A for (A, ν) in zip(table.A, table.NU) if
            even_odd_factors(model, A + ν) == (1.0, 1.0)
        ]
        @test first(unity) == first(masses) && last(unity) == last(masses)
    end

    # ΔZ falls to zero at symmetry where the steep branch reaches ΔZ_max above A_F/2. For CF252S
    # it does not: the steep branch runs down to A_F/2, as Table IV does, and the zero of point X
    # holds at A_F/2 alone.
    U235T = wahl_1988_model("U235T")
    @test charge_polarization(U235T, 118) ≈ 0.0 atol = 1e-12
    @test U235T.A_m > 118
    CF252S = wahl_1988_model("CF252S")
    @test CF252S.A_m < 126
    @test charge_polarization(CF252S, 126) == 0.0
    @test 0.4 < charge_polarization(CF252S, 126.25) < 0.5
    @test charge_dispersion(CF252S, 126) == CF252S.σZ50

    # A fragment of mass A₀/2 is its own complement, so its charge distribution is symmetric
    # about Z₀/2 for every reaction of the evaluation.
    for label in ("U235T", "U233T", "PU239T", "CF252S")
        model = wahl_1988_model(label)
        A = model.A_F ÷ 2
        yields = Dict(zip(fragment_charge_yields(model, A)...))
        @test all(y ≈ get(yields, model.Z_F - Z, 0.0) for (Z, y) in yields)
    end

    # Complementarity, eq. (7d): the light side mirrors the heavy about Z_F/2.
    for A′ in (95.3, 104.0, 110.6)
        @test most_probable_charge(U235T, A′) ≈ 92 - most_probable_charge(U235T, 236 - A′) rtol =
            RTOL
    end
end

@testset "Wahl (1988) against its evaluated experimental data" begin
    # The model is a fit to the evaluated yields whose moments the EXPT columns carry. Where
    # both exist, the model's mean charge sits within a few hundredths of the data on average,
    # and its dispersion within a few hundredths too; these bound the model's representation of
    # the data, not the implementation.
    for (label, table) in WAHL_1988_TABLES
        model = wahl_1988_model(label)
        rows = [i for i in eachindex(table.A) if !isnan(table.AVEZ_expt[i])]
        @test length(rows) >= 9
        Z̄_diff = [
            abs(
                charge_moments(
                    fractional_independent_yields(model, table.A[i], table.NU[i])...,
                )[1] - table.AVEZ_expt[i],
            ) for i in rows
        ]
        @test sum(Z̄_diff) / length(Z̄_diff) < 0.1
    end
end

@testset "fragment-level charge distributions" begin
    model = wahl_1988_model("U235T")
    charges, p = fragment_charge_yields(model, 140)
    @test sum(p) ≈ 1 atol = 1e-12
    @test all(>=(0), p)
    # The proton factor modulates by Z parity: even Z enhanced against the bare lattice Gaussian.
    bare_model = Wahl1988(
        model.label,
        model.channels,
        model.Z_F,
        model.A_F,
        model.ΔZ140,
        model.ΔZSL,
        model.σZ,
        model.σZ50,
        1.0,
        1.0,
        model.ΔA_Z,
        model.ΔZmax,
    )
    _, bare = fragment_charge_yields(bare_model, 140)
    ratio = p ./ bare
    @test all(
        ratio[i] > 1 for i in eachindex(charges) if iseven(charges[i]) && bare[i] > 1e-6
    )
    @test all(
        ratio[i] < 1 for i in eachindex(charges) if isodd(charges[i]) && bare[i] > 1e-6
    )

    # The neutron factor enters only on request, and then on the fragment's own parity.
    paired = Wahl1988("U235T"; neutron_pairing = true)
    @test fragment_charge_yields(paired, 141)[2] != fragment_charge_yields(model, 141)[2]
    @test !model.neutron_pairing

    # The same model through the ChargeModel interface, for the domain.
    system = neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth")
    Zₚ = fragment_most_probable_charge(model, system, 140)
    @test Zₚ ≈ most_probable_charge(model, 140) rtol = RTOL
    @test fragment_charge_probability(model, system, Nuclide(54, 140)) ≈
          p[findfirst(==(54), charges)] rtol = RTOL
    @test fragment_charge_probability(model, system, Nuclide(20, 140)) == 0.0
    @test_throws ArgumentError fragment_charge_probability(model, CF252, Nuclide(54, 140))

    domain = fragmentation_domain(system, model, 118:160)
    @test length(domain) == 43 * 5
    # A five-charge window keeps nearly all of p(Z|A) in the peaks, and visibly less at the
    # Z = 50 crossing, where the distribution is narrow but off-centre from the window.
    kept(A) = sum(e.probability for e in domain.entries if e.heavy.A == A)
    @test kept(140) > 0.999
    @test kept(130) > 0.99

    effective = effective_charge_distribution(model, 118:160)
    @test effective isa ChargeDistribution
    Z̄, RMS = charge_moments(fragment_charge_yields(model, 140)...)
    @test charge_polarization(effective, 140) ≈ Z̄ - 140 * 92 / 236 rtol = RTOL
    @test charge_dispersion(effective, 140) ≈ RMS rtol = RTOL
    @test_throws ArgumentError effective_charge_distribution(model, [100])
end

@testset "the charge model resolves in layers" begin
    table = mass_table()
    thermal = neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth")
    resonance = neutron_induced_fission(Nuclide(92, 235), 5.8013e-4, "nres")
    @test charge_model(table, thermal) isa Wahl1988
    @test charge_model(table, resonance).label == "U235T"
    @test charge_model(table, CF252).label == "CF252S"

    # The nucleus alone is not enough: 240-Pu(sf) reaches 239-Pu(nth,f)'s compound nucleus
    # 6.5 MeV lower, and a fast channel is not the thermal fit.
    @test wahl_1988(spontaneous_fission(Nuclide(94, 240))) === nothing
    @test charge_model(table, spontaneous_fission(Nuclide(94, 240))) isa WahlSystematics
    @test wahl_1988(neutron_induced_fission(Nuclide(92, 235), 2.0, "nfast")) === nothing
    @test charge_model(table, neutron_induced_fission(Nuclide(90, 232), 1.5, "nfast")) isa
          WahlSystematics

    # Outside the systematics, the conventional means.
    outside = charge_model(table, spontaneous_fission(Nuclide(100, 256)))
    @test outside isa ChargeDistribution
    @test outside.source == "mean values"

    # The option reaches the model.
    @test charge_model(table, thermal; neutron_pairing = true).neutron_pairing
    @test_throws ArgumentError Wahl1988("U238F")
end
