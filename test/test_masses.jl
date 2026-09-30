@testset "mass table" begin
    table = mass_table()
    @test length(table) > 3000
    @test haskey(table, NEUTRON)

    # Converted from keV on load, so MeV is the only energy unit inside the package. The value
    # is the evaluation's 8071.31806 keV after the single-precision pass the shipped table went
    # through; data/README.md records the pass.
    neutron = mass_excess(table, NEUTRON)
    @test value(neutron) ≈ 8.07131787109 atol = 1e-9
    @test uncertainty(neutron) ≈ 4.4e-7 atol = 1e-9

    @test mass_excess(table, Nuclide(60, 300)) === nothing
end

@testset "separation energies" begin
    table = mass_table()

    # Sₙ of an actinide sits around 5–7 MeV; a sign or unit slip would leave this range.
    Sₙ = value(neutron_separation_energy(table, Nuclide(98, 252)))
    @test 4.0 < Sₙ < 8.0

    @test neutron_separation_energy(table, Nuclide(60, 300)) === nothing
    @test_throws ArgumentError separation_energy(table, NEUTRON, Nuclide(50, 120))
end

@testset "Q values against the archived run" begin
    table = mass_table()
    reference = read_reference("Q_vs_A_H_Z_H.dat", [:A_H, :Z_H, :Q], 2)
    @test length(reference) == 245

    computed = Dict(
        (A_heavy, Z_heavy) => q_value(table, CF252, Nuclide(Z_heavy, A_heavy)) for
        (A_heavy, Z_heavy) in keys(reference)
    )
    @test count(isnothing, values(computed)) == 0
    # The reference was written at full Float64 precision, so agreement is to round-off.
    @test maximum(abs(value(computed[k]) - reference[k]) for k in keys(reference)) < 1e-9

    # Q for an actinide split is ≈ 200 MeV; every entry must be positive and in range.
    Q_min, Q_max = extrema(value(computed[k]) for k in keys(reference))
    @test Q_min > 150.0
    @test Q_max < 260.0
end

@testset "Q value uncertainties" begin
    table = mass_table()
    heavy = Nuclide(52, 134)
    Q = q_value(table, CF252, heavy)

    # The uncertainty comes from the three mass excesses and nothing else, so it is small
    # against Q itself but non-zero.
    @test uncertainty(Q) > 0
    @test uncertainty(Q) < 0.1

    # Complementarity: naming either fragment of the pair describes the same split.
    light = complementary_fragment(CF252, heavy)
    @test light == Nuclide(46, 118)
    @test value(q_value(table, CF252, light)) ≈ value(Q) rtol = RTOL

    # A fragment heavier or more charged than the compound nucleus is a caller error, not an
    # absent measurement, and is named as such rather than returned as `nothing`.
    @test_throws ArgumentError q_value(table, CF252, Nuclide(60, 300))
    @test_throws ArgumentError complementary_fragment(CF252, Nuclide(98, 252))

    # A physically shaped fragment whose mass excess is simply not tabulated yields `nothing`.
    absent = first(
        nuclide for nuclide in (Nuclide(Z, A) for Z in 20:40, A in 90:110) if
        !haskey(table, nuclide)
    )
    @test q_value(table, CF252, absent) === nothing
end
