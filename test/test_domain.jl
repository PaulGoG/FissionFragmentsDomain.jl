@testset "fragmentation domain" begin
    domain = cf252_domain()

    # 49 heavy masses × 5 charge numbers.
    @test length(domain) == 245
    @test domain.heavy_masses == 126:174
    @test domain.charges_per_mass == 5
    @test mass_range(domain) == 78:174

    # Sorted as whole rows, by heavy mass then heavy charge. Sorting the mass and charge columns
    # independently and looking values up by the pair can report a mass against a charge that
    # never occurred with it.
    keys = [(entry.heavy.A, entry.heavy.Z) for entry in domain]
    @test issorted(keys)
    @test allunique(keys)

    # Mass and charge are conserved by every split, and every entry carries a probability.
    @test all(entry.heavy.A + entry.light.A == 252 for entry in domain)
    @test all(entry.heavy.Z + entry.light.Z == 98 for entry in domain)
    @test minimum(entry.probability for entry in domain) > 0

    # The probability is a property of the pair. For an even-even compound nucleus the light
    # fragment has the parities of the heavy one and the complementary Zₚ, so both carry the
    # same p(Z | A). Away from A₀/2 only: at the symmetric mass the light fragment is evaluated
    # as heavy, and the 1988 CF252S model keeps ΔZ ≈ 0.49 there.
    model = cf252_charge_distribution()
    for entry in (domain[6], domain[120])
        @test entry.probability ≈ fragment_charge_probability(model, CF252, entry.heavy) rtol =
            RTOL
        @test entry.probability ≈ fragment_charge_probability(model, CF252, entry.light) rtol =
            RTOL
    end
end

@testset "symmetric splits" begin
    domain = cf252_domain()
    symmetric = filter(is_symmetric, domain.entries)

    # A_H = 126 is exactly symmetric for 252Cf, so five entries qualify — one per charge.
    @test length(symmetric) == 5
    for entry in symmetric
        @test entry.heavy.A == entry.light.A == 126
    end

    # The fragment view is the nuclide set, each once — what per-nuclide quantities are memoised
    # over. No probability: a nuclide's weight depends on the question being asked.
    listed = fragments(domain)
    @test allunique(listed)
    @test issorted([(nuclide.A, nuclide.Z) for nuclide in listed])
    @test length(listed) == 485
    @test Nuclide(50, 126) in listed
    @test eltype(listed) == Nuclide

    # It agrees with the charge-probability table, the other per-nuclide view.
    @test Set(listed) ==
          Set(keys(fragment_charge_probabilities(domain, cf252_charge_distribution())))
end

@testset "distinct splits" begin
    domain = cf252_domain()
    splits = distinct_splits(domain)

    # Away from mass symmetry a split is labelled once, so only A₀/2 loses entries. The five
    # charge labels at A = 126 describe three physical splits: {47,51}, {48,50} and {49,49}.
    @test length(splits) == length(domain) - 2
    symmetric = [entry for entry in splits if entry.heavy.A == 126]
    @test length(symmetric) == 3
    @test sort([entry.heavy.Z for entry in symmetric]) == [49, 50, 51]

    # Every physical pair appears exactly once, in one labelling or the other.
    pairs = Set(Set([entry.heavy.Z, entry.light.Z]) for entry in symmetric)
    @test length(pairs) == 3

    # A split whose mirror charge falls outside the retained window has no second labelling and
    # must survive whichever side names it. With the mean distribution Zₚ(126) = 48.5, so the
    # window 46:50 is not symmetric about Z₀/2 = 49 — exactly that case.
    offset = fragmentation_domain(CF252, mean_charge_distribution(), 126:174)
    offset_splits = distinct_splits(offset)
    @test length(offset_splits) == length(offset) - 1
    @test sort([e.heavy.Z for e in offset_splits if e.heavy.A == 126]) == [46, 47, 49, 50]
end

@testset "domain validation" begin
    distribution = mean_charge_distribution()

    # The range is validated against the compound nucleus, which is what splits. The
    # prototype validates against the target, which agrees only when the bound sits at A₀/2.
    pu = neutron_induced_fission(Nuclide(94, 239), 2.53e-8, "nth")
    @test symmetric_mass(pu) == 120.0
    @test_throws ArgumentError fragmentation_domain(pu, distribution, 119:160)
    @test length(fragmentation_domain(pu, distribution, 120:160)) == 41 * 5

    @test_throws ArgumentError fragmentation_domain(CF252, distribution, 125:174)
    @test_throws ArgumentError fragmentation_domain(CF252, distribution, 126:252)
    @test_throws ArgumentError fragmentation_domain(
        CF252,
        distribution,
        126:174;
        charges_per_mass = 4,
    )

    # A degenerate single-mass range. The prototype rejects A_H_min >= A_H_max outright.
    @test_throws ArgumentError fragmentation_domain(CF252, distribution, 130:130)

    # Enough charge numbers per mass and the window walks off the end of the charge range,
    # asking for a fragment with no protons or with all of them.
    @test_throws ArgumentError fragmentation_domain(
        CF252,
        distribution,
        126:174;
        charges_per_mass = 199,
    )
end

@testset "split weights" begin
    # The fragment yield counts both fragments of a symmetric split on the heavy side, so a sum
    # over heavy fragments would give that split twice the weight of an asymmetric one.
    @test split_weight(fragmentation_at(126, 49)) == 0.5
    @test split_weight(fragmentation_at(126, 50)) == 0.5
    @test split_weight(fragmentation_at(127, 50)) == 1.0

    # Only the split into two copies of one nuclide is self-complementary.
    @test is_self_complementary(fragmentation_at(126, 49))
    @test !is_self_complementary(fragmentation_at(126, 50))
    @test !is_self_complementary(fragmentation_at(127, 50))
end

@testset "the polarization vanishes at symmetry on request" begin
    # With ΔZ = −0.5 at A₀/2 = 126, Zₚ = 48.5 and the window 46:50 is not its own mirror about
    # Z₀/2 = 49; with the polarization taken as zero there the window is 47:51 and is.
    given = fragmentation_domain(CF252, mean_charge_distribution(), 126:174)
    @test charges(given, 126) == 46:50
    @test symmetric_charge_set_is_invariant(given) == false

    zeroed = fragmentation_domain(
        CF252,
        mean_charge_distribution(),
        126:174;
        zero_polarization_at_symmetry = true,
    )
    @test charges(zeroed, 126) == 47:51
    @test symmetric_charge_set_is_invariant(zeroed) == true

    # Only the symmetric mass is affected.
    @test [e for e in zeroed.entries if e.heavy.A > 126] == [e for e in given.entries if e.heavy.A > 126]

    # Light masses are reached through the complements, and a mass outside the domain has none.
    @test charges(given, 78) == sort([98 - Z for Z in charges(given, 174)])
    @test isempty(charges(given, 60))

    # No symmetric split for an odd compound mass.
    odd = fragmentation_domain(
        neutron_induced_fission(Nuclide(92, 232), 2.53e-8, "nth"),
        mean_charge_distribution(),
        117:150,
    )
    @test symmetric_charge_set_is_invariant(odd) === nothing
end
