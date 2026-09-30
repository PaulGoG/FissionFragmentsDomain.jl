@testset "nuclides" begin
    @test Nuclide(98, 252).Z == 98
    @test neutron_number(Nuclide(98, 252)) == 154
    @test NEUTRON == Nuclide(0, 1)
    @test element_symbol(98) == "Cf"
    @test element_symbol(0) == "n"
    @test element_symbol(118) == "Og"
    @test_throws ArgumentError element_symbol(119)
    @test_throws ArgumentError element_symbol(-1)

    @test_throws ArgumentError Nuclide(-1, 10)
    # A nuclide cannot hold fewer nucleons than protons.
    @test_throws ArgumentError Nuclide(50, 40)
end

@testset "fissioning systems" begin
    cf = spontaneous_fission(Nuclide(98, 252))
    @test is_spontaneous(cf)
    @test cf.channel == "sf"
    @test reaction(cf) == "0,f"
    @test cf.compound == cf.target
    @test iszero(cf.incident_energy)
    @test symmetric_mass(cf) == 126.0

    # The compound nucleus carries one more mass unit than the target, and both are kept:
    # the domain is built on the compound, a fission case is named by the target.
    pu = neutron_induced_fission(Nuclide(94, 239), 2.53e-8, "nth")
    @test pu.compound == Nuclide(94, 240)
    @test pu.target == Nuclide(94, 239)
    @test !is_spontaneous(pu)
    @test reaction(pu) == "n,f"
    @test symmetric_mass(pu) == 120.0

    @test_throws ArgumentError FissioningSystem(
        Nuclide(98, 252),
        Nuclide(98, 252),
        1.0,
        "sf",
    )
    @test_throws ArgumentError FissioningSystem(
        Nuclide(98, 253),
        Nuclide(98, 252),
        0.0,
        "sf",
    )
    @test_throws ArgumentError neutron_induced_fission(Nuclide(94, 239), -1.0, "nth")
    # A channel outside the vocabulary is refused rather than carried into a path.
    @test_throws ArgumentError neutron_induced_fission(Nuclide(94, 239), 0.0, "thermal")
    @test_throws ArgumentError neutron_induced_fission(Nuclide(94, 239), 0.0, "sf")

    @test occursin("252Cf(sf)", sprint(show, cf))
    @test occursin("239Pu(nth,f)", sprint(show, pu))
end

@testset "system label and notation" begin
    # The token names a directory and a configuration file; the notation is for a figure. One
    # name may not mean both, so they are separate functions.
    cf = spontaneous_fission(Nuclide(98, 252))
    @test system_label(cf) == "Cf252_sf"
    @test system_notation(cf) == "²⁵²Cf(sf)"

    thermal = neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth")
    resonance = neutron_induced_fission(Nuclide(92, 235), 5.8013e-4, "nres")
    # The channel is exactly what separates two runs of one target that `n,f` alone cannot.
    @test system_label(thermal) == "U235_nth"
    @test system_label(resonance) == "U235_nres"
    @test reaction(thermal) == reaction(resonance)
    @test system_notation(thermal) == "²³⁵U(nth,f)"
    @test system_notation(resonance) == "²³⁵U(nres,f)"

    @test system_label(neutron_induced_fission(Nuclide(92, 233), 2.53e-8, "nth")) ==
          "U233_nth"
    @test system_label(neutron_induced_fission(Nuclide(94, 239), 2.53e-8, "nth")) ==
          "Pu239_nth"
end

@testset "the system record" begin
    record = system_record(neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth"))
    @test record["label"] == "U235_nth"
    @test record["reaction"] == "n,f"
    @test (record["compound_A"], record["compound_Z"]) == (236, 92)
    @test record["incident_energy_MeV"] == 2.53e-8
    @test length(record) == 9
end
