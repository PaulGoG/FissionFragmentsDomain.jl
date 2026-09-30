using Test
using TOML: TOML
using Measurements: Measurement, measurement, value, uncertainty
using SpecialFunctions: erf
using FissionFragmentsDomain

const REFERENCE = joinpath(@__DIR__, "references", "Cf252_sf")

include("physics.jl")

@testset "FissionFragmentsDomain" begin
    include("test_nuclides.jl")
    include("test_tables.jl")
    include("test_masses.jl")
    include("test_charge.jl")
    include("test_wahl.jl")
    include("test_wahl1988.jl")
    include("test_domain.jl")
    include("test_leveldensity.jl")
    include("test_energetics.jl")
    include("test_temperature_ratio.jl")
    include("test_yields.jl")
    include("test_averaging.jl")
    include("test_curves.jl")
    include("test_manifest.jl")

    @testset "quality" begin
        using Aqua
        Aqua.test_all(FissionFragmentsDomain)

        using ExplicitImports
        @test check_no_implicit_imports(FissionFragmentsDomain) === nothing
        @test check_no_stale_explicit_imports(FissionFragmentsDomain) === nothing

        using JET
        using JET: AnyFrameModule
        using Measurements: Measurements
        # Measurements.jl keeps per-term derivatives in a dictionary whose iteration JET
        # cannot follow; the reports raised there are about that dependency's internals, not
        # about this package, so they are excluded by frame rather than silenced wholesale.
        JET.test_package(
            FissionFragmentsDomain;
            ignored_modules = (AnyFrameModule(Measurements),),
        )
    end
end
