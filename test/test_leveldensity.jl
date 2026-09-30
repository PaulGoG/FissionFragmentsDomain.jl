guarded(SHELL_CORRECTIONS_AVAILABLE, "shell correction file layout") do
    @testset "shell correction file layout" begin
        corrections = read_shell_correction_table(SHELL_CORRECTION_FILE)

        # The reader takes columns by position, so this guards the shipped file rather than the
        # reader: the neutron correction comes first, `n S_N S_Z`. The prototype declares the
        # two the other way round, and because S(Z) is looked up at the proton number while S(N)
        # is looked up at the neutron number, the swap is not absorbed by their sum.
        header = split(strip(first(readlines(SHELL_CORRECTION_FILE))))
        @test header[2] == "S_N"
        @test header[3] == "S_Z"
        @test corrections.S_N[11] ≈ 6.80 rtol = RTOL
        @test corrections.S_Z[11] ≈ -2.91 rtol = RTOL
    end
end

@testset "back-shifted Fermi gas against the archived run" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)
    reference = read_delimited_table(
        joinpath(REFERENCE, "a_vs_A_H_Z_H.dat"),
        TableSpec([:A_H, :Z_H, :a_L, :a_H]; skip = 1),
    )
    heavy_masses = column(reference, :A_H, Int)
    heavy_charges = column(reference, :Z_H, Int)
    a_light = column(reference, :a_L, Float64)
    a_heavy = column(reference, :a_H, Float64)
    @test length(reference) == 244

    heavy = [Nuclide(Z, A) for (A, Z) in zip(heavy_masses, heavy_charges)]
    computed_heavy = [level_density_parameter(model, nuclide) for nuclide in heavy]
    computed_light = [
        level_density_parameter(model, complementary_fragment(CF252, nuclide)) for
        nuclide in heavy
    ]
    @test count(isnothing, computed_heavy) == 0
    @test count(isnothing, computed_light) == 0
    # The reference is rounded to eight decimals on write.
    @test maximum(abs(a - r) for (a, r) in zip(computed_heavy, a_heavy)) < 1e-8
    @test maximum(abs(a - r) for (a, r) in zip(computed_light, a_light)) < 1e-8
end

@testset "BSFG pinned to von Egidy & Bucurescu" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)

    # Eq. (19) of the 2009 paper, coefficients and all. Reconstructing a from S′ by hand must
    # reproduce what the implementation returns, which pins the additive form: the variant
    # p₁(1 + p₂ S′)A^p₃ would differ by a factor p₁ on the S′ term.
    for nuclide in (Nuclide(52, 134), Nuclide(42, 104), Nuclide(50, 132))
        S = shell_correction(model, nuclide)
        @test S !== nothing
        @test level_density_parameter(model, nuclide) ≈
              (0.199 + 0.0096 * S) * nuclide.A^0.869 rtol = RTOL
        @test !isapprox(
            level_density_parameter(model, nuclide),
            0.199 * (1 + 0.0096 * S) * nuclide.A^0.869;
            rtol = RTOL,
        )
    end

    # Eq. (12) of the 2009 paper: the pairing term enters as 0.5 P_a′ with
    # P_a′ = ½[M(A+2,Z+1) − 2M(A,Z) + M(A−2,Z−1)], i.e. a quarter of the second difference,
    # and with no (−1)^Z alternation — that factor distinguishes P_a′ from the P_a of the
    # 2005 paper and of Audi's tables, and it is P_a′ that enters eq. (19).
    nuclide = Nuclide(52, 134)
    Δ = value(mass_excess(table, nuclide))
    Δ_up = value(mass_excess(table, Nuclide(53, 136)))
    Δ_down = value(mass_excess(table, Nuclide(51, 132)))
    pairing = (Δ_up - 2Δ + Δ_down) / 4
    δW₀ =
        liquid_drop_energy(model, nuclide) - (
            nuclide.Z * value(mass_excess(table, Nuclide(1, 1))) +
            neutron_number(nuclide) * value(mass_excess(table, NEUTRON)) - Δ
        )
    @test shell_correction(model, nuclide) ≈ δW₀ + pairing rtol = RTOL

    # S = M_exp − M_LD (2005, eq. 7), so a nucleus more bound than the liquid drop has S < 0.
    # Shell closures are exactly those nuclei, and the suppressed a is what puts the sawtooth
    # into E*(A) and hence into ν(A).
    @test shell_correction(model, Nuclide(50, 132)) < -8.0      # doubly magic
    @test shell_correction(model, Nuclide(42, 104)) > 0.0       # mid-shell
    @test abs(shell_correction(model, Nuclide(50, 120))) < 1.0  # magic Z, mid-shell N

    # a/A ≈ 1/8 away from closures, and well below it at one.
    @test level_density_parameter(model, Nuclide(42, 104)) / 104 > 0.11
    @test level_density_parameter(model, Nuclide(50, 132)) / 132 < 0.06
end

@testset "level density parameter behaviour" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)

    # Roughly A/8 MeV⁻¹ across the fragment range, the familiar magnitude for these nuclei.
    for nuclide in (Nuclide(52, 134), Nuclide(46, 118), Nuclide(38, 96))
        a = level_density_parameter(model, nuclide)
        @test a !== nothing
        @test 0.05 * nuclide.A < a < 0.2 * nuclide.A
    end

    # The N = 82, Z = 50 closure suppresses the parameter relative to its neighbours, which is
    # the origin of the sawtooth in E*(A) and hence in ν(A).
    magic = level_density_parameter(model, Nuclide(50, 132))
    neighbour = level_density_parameter(model, Nuclide(50, 138))
    @test magic < neighbour

    # Unavailable rather than an error where the neighbouring masses are not tabulated.
    @test level_density_parameter(model, Nuclide(1, 1)) === nothing

    if SHELL_CORRECTIONS_AVAILABLE
        gilbert = GilbertCameron(read_shell_correction_table(SHELL_CORRECTION_FILE))
        a_gc = level_density_parameter(gilbert, Nuclide(52, 134))
        @test a_gc !== nothing
        @test 0.05 * 134 < a_gc < 0.2 * 134

        # Outside the tabulated nucleon numbers the systematics simply has nothing to say.
        @test level_density_parameter(gilbert, Nuclide(2, 4)) === nothing
    end
end

@testset "the Gilbert-Cameron deformed region" begin
    # The paper's own definition, p. 1459: two rectangles in (Z, N), bounds inclusive.
    @test is_deformed_gilbert_cameron(Nuclide(58, 148))    # Z=58, N=90 — heavy fragment peak
    @test is_deformed_gilbert_cameron(Nuclide(54, 140))    # on the lower corner, Z=54, N=86
    @test is_deformed_gilbert_cameron(Nuclide(78, 200))    # upper corner of the first region
    @test !is_deformed_gilbert_cameron(Nuclide(52, 134))   # Z=52 below it, N=82 magic
    @test !is_deformed_gilbert_cameron(Nuclide(40, 100))   # light fragment, undeformed
    @test !is_deformed_gilbert_cameron(Nuclide(54, 136))   # Z in range, N=82 is not
    @test is_deformed_gilbert_cameron(Nuclide(94, 239))    # second region
    @test is_deformed_gilbert_cameron(Nuclide(15, 31))     # light nuclei, p. 1464
    @test !is_deformed_gilbert_cameron(Nuclide(25, 55))    # 20 ≤ Z ≤ 29 follows line I
end

guarded(SHELL_CORRECTIONS_AVAILABLE, "the two Gilbert-Cameron branches") do
    @testset "the two Gilbert-Cameron branches" begin
        gilbert = GilbertCameron(read_shell_correction_table(SHELL_CORRECTION_FILE))

        # The branches differ by the intercept alone, 0.142 − 0.120 = 0.022 per nucleon, so the
        # deformed one is lower by exactly 0.022 A. Checked on a fragment inside the region.
        deformed = Nuclide(58, 148)
        @test is_deformed_gilbert_cameron(deformed)
        a = level_density_parameter(gilbert, deformed)
        corrections = read_shell_correction_table(SHELL_CORRECTION_FILE)
        S = corrections.S_Z[58] + corrections.S_N[90]
        undeformed_value = 148 * (0.00917 * S + 0.142)
        @test undeformed_value - a ≈ 0.022 * 148 rtol = RTOL

        # Roughly a fifth: the size of the correction on the heavy fragment peak.
        @test 0.15 < (undeformed_value - a) / a < 0.35

        # Selecting eq. (20) alone reverts to the single-branch behaviour, which is what the
        # published results of this method used and what the companion package implements.
        undeformed = GilbertCameron(
            read_shell_correction_table(SHELL_CORRECTION_FILE);
            deformed_branch = false,
        )
        @test level_density_parameter(undeformed, deformed) ≈ undeformed_value rtol = RTOL
        @test level_density_parameter(undeformed, Nuclide(40, 100)) ≈
              level_density_parameter(gilbert, Nuclide(40, 100)) rtol = RTOL

        # Z ≥ 99 has no tabulated S(Z). The file pads it with zeros; the reader must not store
        # them, or the systematics silently answers where the paper says nothing.
        # Both have a tabulated S(N); they differ only in whether S(Z) exists.
        @test level_density_parameter(gilbert, Nuclide(99, 240)) === nothing
        @test level_density_parameter(gilbert, Nuclide(98, 240)) !== nothing
    end
end
