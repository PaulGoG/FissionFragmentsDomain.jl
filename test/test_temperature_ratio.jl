@testset "temperature ratio and the excitation-energy partition" begin
    # The forward relation and its inverse close on each other.
    for r_ν in (0.30, 0.45, 0.62), R_a in (0.85, 1.0, 1.2)
        R_T = temperature_ratio(r_ν, R_a)
        @test heavy_excitation_fraction(R_T, R_a) ≈ r_ν rtol = RTOL
        # E*_L/E*_H = R_a R_T² = (1 − r_ν)/r_ν.
        @test R_a * R_T^2 ≈ (1 - r_ν) / r_ν rtol = RTOL
    end
    @test round(temperature_ratio(0.40, 1.05); digits = 4) == 1.1952
    @test heavy_excitation_fraction(1.0, 1.0) == 0.5

    # The slope is the derivative of the relation.
    r_ν, R_a, h = 0.42, 1.07, 1e-6
    R_T = temperature_ratio(r_ν, R_a)
    numerical = (temperature_ratio(r_ν + h, R_a) - temperature_ratio(r_ν - h, R_a)) / (2h)
    @test temperature_ratio_slope(R_T, R_a, r_ν) ≈ numerical rtol = 1e-7

    @test_throws DomainError temperature_ratio(0.0, 1.0)
    @test_throws DomainError temperature_ratio(1.0, 1.0)
    @test_throws DomainError temperature_ratio(0.5, 0.0)
    @test_throws DomainError heavy_excitation_fraction(0.0, 1.0)
    @test_throws DomainError heavy_excitation_fraction(1.0, -1.0)
end

@testset "level density ratio over the charge distribution" begin
    model = BackShiftedFermiGas(mass_table())
    domain = fragmentation_domain(
        CF252,
        mean_charge_distribution(),
        CF252_HEAVY_MASSES;
        zero_polarization_at_symmetry = true,
    )
    of_means = level_density_ratio(RatioOfMeans(), model, domain)
    of_ratios = level_density_ratio(MeanOfRatios(), model, domain)
    @test keys(of_means) == keys(of_ratios) == Set(CF252_HEAVY_MASSES)

    # At symmetry, with a mirror-invariant charge set, the two averages of the ratio of means
    # run over the same nuclides: R_a = 1 up to the order of summation. The mean of ratios
    # averages reciprocal pairs and lies above one.
    @test of_means[126] ≈ 1 atol = 4 * eps()
    @test of_ratios[126] > 1 + 1e-4
    @test of_ratios[126] < 1 + 1e-2

    # Away from symmetry the two orders differ at second order in the spread of a over the
    # charge window, which is largest where a changes fastest with Z: at the doubly magic heavy
    # fragment, ¹³²Sn.
    gap = Dict(A => abs(of_ratios[A] / of_means[A] - 1) for A in 127:174)
    @test maximum(values(gap)) < 0.02
    @test argmax(gap) in 130:134

    # A single-charge window reduces both orders to the ratio of that pair.
    narrow = fragmentation_domain(
        CF252,
        mean_charge_distribution(),
        140:141;
        charges_per_mass = 1,
    )
    entry = first(narrow.entries)
    expected =
        level_density_parameter(model, entry.light) /
        level_density_parameter(model, entry.heavy)
    @test level_density_ratio(RatioOfMeans(), model, narrow)[140] ≈ expected rtol = RTOL
    @test level_density_ratio(MeanOfRatios(), model, narrow)[140] ≈ expected rtol = RTOL
end

@testset "charge-resolved inversion" begin
    model = BackShiftedFermiGas(mass_table())
    domain = fragmentation_domain(CF252, cf252_charge_distribution(), CF252_HEAVY_MASSES)
    resolved = ChargeResolved()

    # The inverse closes on the charge-reduced partition to rounding, over the whole domain.
    for A_H in CF252_HEAVY_MASSES, R_T in (0.8, 1.0, 1.45)
        r_ν = heavy_excitation_fraction(resolved, model, domain, A_H, R_T)
        @test temperature_ratio(resolved, model, domain, A_H, r_ν) ≈ R_T rtol = 1e-13
    end

    # The effective ratios do not close on it: the gap is second order in the spread of a_L/a_H
    # and largest at the doubly magic heavy fragment.
    miss(averaging) = maximum(
        abs(
            temperature_ratio(
                averaging,
                model,
                domain,
                A_H,
                heavy_excitation_fraction(resolved, model, domain, A_H, 1.45),
            ) - 1.45,
        ) for A_H in CF252_HEAVY_MASSES
    )
    @test 1e-4 < miss(RatioOfMeans()) < miss(MeanOfRatios()) < 0.02

    # At the symmetric split the charge set is its own mirror, so r_ν = 1/2 gives R_T = 1.
    @test symmetric_charge_set_is_invariant(domain)
    @test temperature_ratio(resolved, model, domain, 126, 0.5) ≈ 1 atol = 4 * eps()

    # A single charge per mass reduces every treatment to the closed form of that pair.
    narrow = fragmentation_domain(
        CF252,
        cf252_charge_distribution(),
        140:141;
        charges_per_mass = 1,
    )
    entry = first(narrow.entries)
    ρ =
        level_density_parameter(model, entry.light) /
        level_density_parameter(model, entry.heavy)
    for averaging in (resolved, RatioOfMeans(), MeanOfRatios())
        @test temperature_ratio(averaging, model, narrow, 140, 0.4) ≈
              temperature_ratio(0.4, ρ) rtol = 1e-14
    end

    # The slope is the derivative of the inverse.
    A_H, r_ν, h = 132, 0.37, 1e-6
    R_T = temperature_ratio(resolved, model, domain, A_H, r_ν)
    numerical =
        (
            temperature_ratio(resolved, model, domain, A_H, r_ν + h) -
            temperature_ratio(resolved, model, domain, A_H, r_ν - h)
        ) / (2h)
    @test temperature_ratio_slope(resolved, model, domain, A_H, R_T) ≈ numerical rtol = 1e-6
    for averaging in (RatioOfMeans(), MeanOfRatios())
        R = temperature_ratio(averaging, model, domain, A_H, r_ν)
        @test temperature_ratio_slope(averaging, model, domain, A_H, R) ≈
              temperature_ratio_slope(
            R,
            level_density_ratio(averaging, model, domain)[A_H],
            r_ν,
        ) rtol = RTOL
    end

    # A mass outside the domain has no fragmentation; an r_ν outside (0, 1) has no root.
    @test temperature_ratio(resolved, model, domain, 200, 0.4) === nothing
    @test heavy_excitation_fraction(resolved, model, domain, 200, 1.0) === nothing
    @test_throws DomainError temperature_ratio(resolved, model, domain, 140, 1.0)
    @test_throws DomainError heavy_excitation_fraction(resolved, model, domain, 140, 0.0)
    @test_throws ArgumentError level_density_ratio(resolved, model, domain)
end
