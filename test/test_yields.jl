@testset "experimental yield input" begin
    mktempdir() do dir
        path = joinpath(dir, "Y.dat")
        write(
            path,
            """
            A TKE Value σ
            140 170.0 1.5 0.1
            140 171.0 2.5 0.2
            141 170.0 3.0 NaN
            141 171.0 NaN NaN
            """,
        )
        measured = read_mass_energy_yield(path)

        # A cell with no finite yield is dropped: an unmeasured (A, TKE) is not a yield of zero.
        @test length(measured) == 3
        @test measured.values[(140, 171.0)] == 2.5
        @test !haskey(measured.values, (141, 171.0))
        @test measured.masses == [140, 141]
        @test measured.energies == [170.0, 171.0]
        @test occursin("2 masses", sprint(show, measured))

        # The uncertainty column is optional, as everywhere else.
        bare = joinpath(dir, "bare.dat")
        write(bare, "A TKE Value\n140 170.0 1.5\n")
        @test length(read_mass_energy_yield(bare)) == 1

        # A repeat carrying the same yield is the symmetric mass written from both sides of
        # the split, which the 235-U distribution does at every kinetic energy. One measurement,
        # counted once.
        mirrored = joinpath(dir, "mirrored.dat")
        write(mirrored, "A TKE Value\n118 170.0 1.5\n118 170.0 1.5\n119 170.0 2.0\n")
        @test length(read_mass_energy_yield(mirrored)) == 2
        @test read_mass_energy_yield(mirrored).values[(118, 170.0)] == 1.5

        # A repeat that disagrees is a real inconsistency in the file.
        repeated = joinpath(dir, "repeated.dat")
        write(repeated, "A TKE Value\n140 170.0 1.5\n140 170.0 2.0\n")
        @test_throws ArgumentError read_mass_energy_yield(repeated)

        empty = joinpath(dir, "empty.dat")
        write(empty, "A TKE Value\n140 170.0 NaN\n")
        @test_throws ArgumentError read_mass_energy_yield(empty)
    end
end

@testset "coverage reports what a sweep reaches" begin
    yields = FragmentYield(
        Dict((140, 54, 170.0) => 3.0, (140, 54, 171.0) => 1.0, (140, 54, 172.0) => 4.0),
        8.0,
    )
    report = coverage(yields, 170.0:2.0:172.0)
    @test report.reached ≈ 7 / 8 rtol = RTOL
    @test report.missed ≈ 1 / 8 rtol = RTOL
    @test report.tabulated == [170.0, 171.0, 172.0]
    @test report.swept == [170.0, 172.0]
end

@testset "yield construction validation" begin
    @test_throws ArgumentError yield_over(
        FragmentYield(Dict((1, 1, 1.0) => 1.0), 1.0),
        :temperature,
    )
end

@testset "mass-dependent kinetic energy dispersion" begin
    mktempdir() do directory
        # A table that reaches only part of the mass range, so the fallback is exercised too.
        path = joinpath(directory, "sigma_TKE_vs_A.dat")
        write(path, "A sigma_TKE\n120 7.0\n121 13.0\n122 9.0\n")
        dispersion = read_kinetic_energy_dispersion(path; default = 10.0)
        @test kinetic_energy_dispersion(dispersion, 120) ≈ 7.0 rtol = RTOL
        @test kinetic_energy_dispersion(dispersion, 121) ≈ 13.0 rtol = RTOL
        # Outside the table it falls back, exactly as a charge distribution does.
        @test kinetic_energy_dispersion(dispersion, 200) ≈ 10.0 rtol = RTOL
        @test dispersion.source == path

        # Columns are taken by position, so the header text is irrelevant.
        renamed = joinpath(directory, "renamed.dat")
        write(renamed, "mass width\n120 7.0\n121 13.0\n")
        @test kinetic_energy_dispersion(read_kinetic_energy_dispersion(renamed), 120) ≈ 7.0 rtol =
            RTOL

        # A non-positive width is refused rather than producing a degenerate Gaussian.
        bad = joinpath(directory, "bad.dat")
        write(bad, "A sigma_TKE\n120 0.0\n")
        @test_throws ArgumentError read_kinetic_energy_dispersion(bad)
        @test_throws ArgumentError uniform_kinetic_energy_dispersion(-1.0)

        # One width everywhere is the no-table case, and it is what a bare number means.
        uniform = uniform_kinetic_energy_dispersion(8.0)
        @test all(
            isapprox(kinetic_energy_dispersion(uniform, A), 8.0; rtol = RTOL) for
            A in 100:160
        )
    end
end

@testset "the dispersion table reaches the reconstruction" begin
    mktempdir() do directory
        # Two marginals that overlap on three masses, reconstructed on a coarse grid.
        write(joinpath(directory, "Y.dat"), "A Y\n118 10.0\n119 20.0\n120 30.0\n")
        write(joinpath(directory, "TKE.dat"), "A TKE\n118 170.0\n119 170.0\n120 170.0\n")
        grid = 150.0:5.0:190.0

        narrow = factorized_yield(
            joinpath(directory, "Y.dat"),
            joinpath(directory, "TKE.dat"),
            grid,
            236;
            dispersion = 5.0,
        )
        wide = factorized_yield(
            joinpath(directory, "Y.dat"),
            joinpath(directory, "TKE.dat"),
            grid,
            236;
            dispersion = 15.0,
        )
        # A narrower width concentrates each mass row on the mean kinetic energy.
        peak(y, A) = y.values[(A, 170.0)] / sum(y.values[(A, T)] for T in grid)
        @test peak(narrow, 119) > peak(wide, 119)

        # A table giving one mass the narrow width and another the wide one reproduces each.
        write(joinpath(directory, "sigma.dat"), "A sigma_TKE\n118 5.0\n120 15.0\n")
        mixed = factorized_yield(
            joinpath(directory, "Y.dat"),
            joinpath(directory, "TKE.dat"),
            grid,
            236;
            dispersion = read_kinetic_energy_dispersion(
                joinpath(directory, "sigma.dat");
                default = 5.0,
            ),
        )
        @test peak(mixed, 118) ≈ peak(narrow, 118) rtol = RTOL
        @test peak(mixed, 120) ≈ peak(wide, 120) rtol = RTOL
        # 119 is absent from the table and takes the fallback, which here is the narrow
        # width.
        @test peak(mixed, 119) ≈ peak(narrow, 119) rtol = RTOL
    end
end

@testset "one-dimensional mass yield" begin
    mktempdir() do dir
        path = joinpath(dir, "Y_vs_A.dat")
        write(path, "A Y Y_uncertainty\n140 6.1 0.2\n132 4.8 0\n150 1.0 -1\n")
        yields = read_mass_yield(path)
        @test yields.A == [132, 140, 150]
        @test yields.Y == [4.8, 6.1, 1.0]
        @test isequal(yields.σY, [missing, 0.2, missing])
        @test yields.label == "Y_vs_A"
        @test length(yields) == 3
        @test isequal(mass_yield(yields, 132), (4.8, missing))
        @test mass_yield(yields, 133) === nothing

        bare = joinpath(dir, "bare.dat")
        write(bare, "A Y\n140 6.1\n")
        @test isequal(read_mass_yield(bare; label = "x").σY, [missing])

        negative = joinpath(dir, "negative.dat")
        write(negative, "A Y\n140 -6.1\n")
        @test_throws ArgumentError read_mass_yield(negative)

        repeated = joinpath(dir, "repeated.dat")
        write(repeated, "A Y\n140 6.1\n140 6.1\n")
        @test_throws ArgumentError read_mass_yield(repeated)
    end
end
