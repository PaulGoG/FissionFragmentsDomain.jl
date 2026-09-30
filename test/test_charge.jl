guarded(CHARGE_DISTRIBUTION_AVAILABLE, "charge distribution") do
    @testset "charge distribution" begin
        distribution = cf252_charge_distribution()
        @test charge_polarization(distribution, 126) ≈ 0.41847 rtol = RTOL
        @test charge_dispersion(distribution, 126) ≈ 0.50029 rtol = RTOL

        # A mass outside the tabulated range falls back to the conventional means, which is a
        # documented behaviour rather than a silent zero.
        @test charge_polarization(distribution, 300) == -0.5
        @test charge_dispersion(distribution, 300) == 0.6

        mean = mean_charge_distribution()
        @test charge_polarization(mean, 126) == -0.5
        @test charge_dispersion(mean, 126) == 0.6
    end
end

@testset "a mass column written in floating point still reads" begin
    # The 235-U table writes its mass column as 118.0, 119.0, … because a generator emitted it
    # in floating point. Read strictly as Int the file cannot be loaded at all, which blocked
    # 235-U entirely. Checked against the shipped table where it is present, and against a
    # written fixture regardless, so the behaviour is pinned on a bare clone too.
    if U235_POLARIZATION_AVAILABLE
        distribution = read_charge_distribution(U235_CHARGE_DISTRIBUTION_FILE)
        @test charge_polarization(distribution, 118) ≈ 0.00137468 rtol = RTOL
        @test charge_dispersion(distribution, 118) ≈ 0.62081 rtol = RTOL
    end

    mktempdir() do dir
        floating = joinpath(dir, "floating.dat")
        write(floating, "A ΔZ rms\n118.0 0.00137468 0.62081\n119.0 0.0766041 0.6204323\n")
        distribution = read_charge_distribution(floating)
        @test charge_polarization(distribution, 118) ≈ 0.00137468 rtol = RTOL
        @test charge_dispersion(distribution, 118) ≈ 0.62081 rtol = RTOL
    end

    # A genuinely fractional nucleon number is still an error, and says what it could not read.
    mktempdir() do dir
        bad = joinpath(dir, "fractional.dat")
        write(bad, "A ΔZ rms\n118.5 0.1 0.6\n")
        @test_throws ArgumentError read_charge_distribution(bad)
    end
end

@testset "most probable charge" begin
    # Z_UCD(A) = A Z₀/A₀, exact at the symmetric split.
    @test unchanged_charge_distribution(CF252, 126) ≈ 49.0 rtol = RTOL
    @test most_probable_charge(CF252, 126, 0.0) ≈ 49.0 rtol = RTOL
    @test most_probable_charge(CF252, 126, 0.41847) ≈ 49.41847 rtol = RTOL
end

@testset "isobaric charge distribution" begin
    # A normalized Gaussian: the peak sits at 1/(√(2π) rms) and it is symmetric about Zₚ.
    @test charge_probability(50.0, 50.0, 0.6) ≈ 1 / (sqrt(2π) * 0.6) rtol = RTOL
    @test charge_probability(49.0, 50.0, 0.6) ≈ charge_probability(51.0, 50.0, 0.6) rtol =
        RTOL

    # Summed over a wide enough integer range it approaches unity, but not exactly: sampling a
    # Gaussian of width 0.6 on an integer lattice overshoots by 2exp(−2π²rms²) ≈ 1.6e-3 by
    # Poisson summation. The excess is a property of the discretization, not of the expression.
    total = sum(charge_probability(Z, 50.0, 0.6) for Z in 40:60)
    @test total ≈ 1.0 atol = 3e-3
    @test total > 1.0

    @test_throws ArgumentError charge_probability(50.0, 50.0, 0.0)
    @test_throws ArgumentError charge_probability(50.0, 50.0, -0.6)
end

@testset "charge numbers per mass" begin
    @test charge_numbers(49.41847, 5) == 47:51
    @test charge_numbers(49.6, 1) == 50:50
    @test length(charge_numbers(49.4, 3)) == 3

    # An even count cannot be centred on the rounded most probable charge without biasing
    # the domain to one side.
    @test_throws ArgumentError charge_numbers(49.4, 4)
    @test_throws ArgumentError charge_numbers(49.4, 0)
end
