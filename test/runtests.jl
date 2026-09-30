using FissionFragmentsDomain
using Test
using Aqua
using JET

@testset "FissionFragmentsDomain.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(FissionFragmentsDomain)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(FissionFragmentsDomain; target_defined_modules = true)
    end
    # Write your tests here.
end
