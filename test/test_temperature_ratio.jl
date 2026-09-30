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
