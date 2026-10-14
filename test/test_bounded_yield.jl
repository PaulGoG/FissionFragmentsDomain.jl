# Thermal ²³⁵U, heavy masses 130–134 and their complements, the conventional charge means: the
# limits Q + E*_CN of these splits lie between 187 and 206 MeV, so a mean TKE of 185 MeV with a
# 10 MeV width puts a measurable share of every row beyond them.
const U235 = neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth")

function u235_marginals(directory; mean, width)
    masses = vcat(102:106, 130:134)
    write(
        joinpath(directory, "Y.dat"),
        "A Y\n" * join(("$A $(1.0 + 0.1 * (A % 7))" for A in masses), "\n") * "\n",
    )
    write(
        joinpath(directory, "TKE.dat"),
        "A TKE\n" * join(("$A $(mean + 0.3 * (A - 132))" for A in 130:134), "\n") * "\n",
    )
    write(
        joinpath(directory, "sigma.dat"),
        "A sigma_TKE\n" * join(("$A $width" for A in 130:134), "\n") * "\n",
    )
    return joinpath.(directory, ("Y.dat", "TKE.dat", "sigma.dat"))
end

# Mean, standard deviation and yield of every mass of a fragment yield, charges pooled.
function mass_moments(yield::FragmentYield)
    sums = Dict{Int, NTuple{3, Float64}}()
    for ((A, _, TKE), Y) in yield.values
        S = get(sums, A, (0.0, 0.0, 0.0))
        sums[A] = (S[1] + Y, S[2] + Y * TKE, S[3] + Y * TKE^2)
    end
    return Dict(
        A => (S[2] / S[1], sqrt(S[3] / S[1] - (S[2] / S[1])^2), S[1]) for (A, S) in sums
    )
end

@testset "the physical limit of the kinetic energy" begin
    table = mass_table()
    domain = fragmentation_domain(U235, mean_charge_distribution(), 130:134)
    E_compound = value(compound_nucleus_excitation(table, U235))
    heavy = Nuclide(50, 132)
    @test kinetic_energy_limit(table, U235, heavy) ≈
          value(q_value(table, U235, heavy)) + E_compound rtol = 1e-14
    # The limit belongs to the split, so both partners carry it.
    @test kinetic_energy_limit(table, U235, heavy) ==
          kinetic_energy_limit(table, U235, complementary_fragment(U235, heavy))
    largest = maximum_kinetic_energy_limit(table, U235, domain)
    @test largest.limit ≈
          maximum(kinetic_energy_limit(table, U235, f.heavy) for f in domain)
    @test largest.fragmentation.heavy == heavy
    # No configuration at the limit has excitation energy left.
    @test energetics(table, U235, largest.fragmentation, largest.limit + 1e-9) === nothing
end

@testset "a bounded reconstruction keeps the measured moments" begin
    table = mass_table()
    distribution = mean_charge_distribution()
    domain = fragmentation_domain(U235, distribution, 130:134)
    grid = 140.0:1.0:210.0
    mktempdir() do directory
        Y_path, TKE_path, σ_path = u235_marginals(directory; mean = 185.0, width = 10.0)
        dispersion = read_kinetic_energy_dispersion(σ_path)
        build(moments) = bounded_fragment_yield(
            domain,
            distribution,
            table,
            Y_path,
            TKE_path,
            grid;
            dispersion = dispersion,
            symmetrize = true,
            moments = moments,
        )
        measured = kinetic_energy_moments(
            Y_path,
            TKE_path,
            236;
            dispersion = dispersion,
            symmetrize = true,
        )
        plain = build_fragment_yield(
            domain,
            distribution,
            factorized_yield(
                Y_path,
                TKE_path,
                grid,
                236;
                dispersion = dispersion,
                symmetrize = true,
            ),
        )
        plain_moments = mass_moments(plain)

        recentred, bounds = build(:mean)
        @test sum(recentred) ≈ DEFAULT_YIELD_TOTAL rtol = 1e-12
        # Nothing at or beyond the limit survives.
        @test all(
            TKE < kinetic_energy_limit(table, U235, Nuclide(Z, A)) for
            (A, Z, TKE) in keys(recentred.values)
        )
        # The cut is active: the unbounded rows reach past the limits.
        @test all(bound.removed > 1e-3 for bound in bounds)
        moments = mass_moments(recentred)
        for bound in bounds
            A = bound.A
            target = first(measured[A])
            # The truncated mean is the measured one; the plain construction falls below it once
            # its cells beyond the limit are left out.
            @test moments[A][1] ≈ target atol = 1e-9
            @test bound.mean ≈ target atol = 1e-9
            @test bound.measured_mean == target
            # Re-centred upward with the width kept: the cut removed the upper tail.
            @test bound.parent_mean > target
            @test bound.parent_dispersion ≈ 10.0 rtol = 1e-12
            @test bound.dispersion < 10.0
            # The mass yield is the measured one, whatever the cut removed.
            @test moments[A][3] / plain_moments[A][3] ≈ 1 rtol = 1e-12
        end
        # The charge distribution is kept: Y(A, Z) ∝ Y(A) p(Z, A) at every retained charge.
        charge_yield(y) = Dict(
            (A, Z) => sum(v for ((B, W, _), v) in y.values if B == A && W == Z) for
            (A, Z, _) in keys(y.values)
        )
        bounded_charges = charge_yield(recentred)
        plain_charges = charge_yield(plain)
        @test all(
            isapprox(bounded_charges[key] / plain_charges[key], 1; rtol = 1e-12) for
            key in keys(bounded_charges)
        )

        # A bound without the moments restored lowers every mean the cut reaches.
        cut, cut_bounds = build(:none)
        for bound in cut_bounds
            @test bound.mean < bound.measured_mean
            @test mass_moments(cut)[bound.A][1] ≈ bound.mean atol = 1e-9
        end

        # Restoring the width as well: mean and standard deviation both the measured ones.
        rescaled, rescaled_bounds = build(:mean_and_width)
        for bound in rescaled_bounds
            @test mass_moments(rescaled)[bound.A][1] ≈ first(measured[bound.A]) atol = 1e-9
            @test mass_moments(rescaled)[bound.A][2] ≈ 10.0 atol = 1e-8
            @test bound.parent_dispersion > 10.0
        end

        # A grid that stops short of a limit would cut the distribution where it is physical.
        @test_throws ArgumentError bounded_fragment_yield(
            domain,
            distribution,
            table,
            Y_path,
            TKE_path,
            140.0:1.0:200.0;
            dispersion = dispersion,
            symmetrize = true,
        )
        @test_throws ArgumentError build(:median)
    end
end

@testset "without a cut the bounded reconstruction is the plain Gaussian" begin
    table = mass_table()
    distribution = mean_charge_distribution()
    domain = fragmentation_domain(U235, distribution, 130:134)
    grid = 100.0:1.0:210.0
    mktempdir() do directory
        # A mean 45 MeV below the lowest limit with a 5 MeV width: nothing measurable is cut.
        Y_path, TKE_path, σ_path = u235_marginals(directory; mean = 140.0, width = 5.0)
        dispersion = read_kinetic_energy_dispersion(σ_path)
        plain = build_fragment_yield(
            domain,
            distribution,
            factorized_yield(
                Y_path,
                TKE_path,
                grid,
                236;
                dispersion = dispersion,
                symmetrize = true,
            ),
        )
        for moments in KINETIC_ENERGY_MOMENTS
            bounded, bounds = bounded_fragment_yield(
                domain,
                distribution,
                table,
                Y_path,
                TKE_path,
                grid;
                dispersion = dispersion,
                symmetrize = true,
                moments = moments,
            )
            @test all(bound.removed < 1e-15 for bound in bounds)
            @test all(
                abs(bound.parent_mean - bound.measured_mean) < 1e-9 for bound in bounds
            )
            @test all(bound.parent_dispersion ≈ 5.0 for bound in bounds)
            # Cell by cell the plain Gaussian; the plain cells beyond the limit carry nothing.
            @test all(
                isapprox(bounded.values[key], cell; rtol = 1e-10) for
                (key, cell) in plain.values if haskey(bounded.values, key)
            )
            @test sum(
                (cell for (key, cell) in plain.values if !haskey(bounded.values, key));
                init = 0.0,
            ) < 1e-12
        end
    end
end

@testset "a bounded joint distribution keeps its own row moments" begin
    table = mass_table()
    distribution = mean_charge_distribution()
    domain = fragmentation_domain(U235, distribution, 130:134)
    grid = 160.0:2.0:210.0
    # A measured matrix, skewed, part of each row beyond the limits.
    values = Dict{Tuple{Int, Float64}, Float64}()
    for A in vcat(102:106, 130:134), TKE in grid
        values[(A, TKE)] = (TKE - 155.0) * exp(-(TKE - 185.0)^2 / 200.0)
    end
    joint = MassEnergyYield(
        values,
        sort(unique(first.(keys(values)))),
        collect(grid),
        "fixture",
    )
    row_mean(A) =
        sum(TKE * values[(A, TKE)] for TKE in grid) / sum(values[(A, TKE)] for TKE in grid)
    bounded, bounds = bounded_fragment_yield(domain, distribution, table, joint)
    moments = mass_moments(bounded)
    for bound in bounds
        @test bound.measured_mean ≈ row_mean(bound.A) rtol = 1e-14
        @test moments[bound.A][1] ≈ row_mean(bound.A) atol = 1e-9
        @test isnan(bound.parent_mean)
        @test bound.removed > 0
    end
    # Without the moments restored, each charge keeps the measured row below its limit, in the
    # measured proportions.
    cut, _ = bounded_fragment_yield(domain, distribution, table, joint; moments = :none)
    for (A, Z) in unique((A, Z) for (A, Z, _) in keys(cut.values))
        ratios =
            [Y / values[(A, TKE)] for ((B, W, TKE), Y) in cut.values if B == A && W == Z]
        @test maximum(ratios) ≈ minimum(ratios) rtol = 1e-12
    end
end
