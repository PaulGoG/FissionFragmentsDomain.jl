@testset "Wahl Zₚ systematics" begin
    # 235-U(n_th,f), the reaction the systematics is referenced to: the differences in eq. (17)
    # all vanish except the excitation term, so the peak parameters reduce to their intercepts.
    U236 = WahlSystematics(92, 236, 6.551)
    @test U236.σZ140 ≈ 0.566 rtol = RTOL
    @test U236.ΔZ140 ≈ -0.487 rtol = RTOL
    @test U236.σZ50 ≈ 0.356 rtol = RTOL
    @test U236.SL50 ≈ 0.191 rtol = RTOL
    @test U236.ΔZmax ≈ 0.699 rtol = RTOL

    # Region boundaries, eq. (10) and Fig. 18, p. 31: the wings start at 77, the near-symmetry
    # region spans roughly 106 to 130, and the Z = 50 crossing sits near 125.
    @test U236.B2 ≈ 77.0 rtol = RTOL
    @test 129 < U236.B4 < 132
    @test 110 < U236.Ba < 113
    @test 123 < U236.Bb < 126
    @test U236.B5 ≈ 236 - U236.B2 rtol = RTOL
    @test U236.B6 ≈ 166.0 rtol = RTOL

    # Charge conservation: the distribution is equal and opposite about the symmetric split, so
    # it vanishes there. Nothing in the construction imposes this — it follows from the straight
    # line between Ba and Bb being centred on A_F/2.
    @test wahl_polarization(U236, 118) ≈ 0.0 atol = 1e-12

    # ΔZ is continuous across every boundary the heavy side crosses: the legs are constructed to
    # meet, and a jump in the centroid of the charge distribution would be unphysical.
    for boundary in (U236.Bb, U236.B4, U236.B5, U236.B6)
        below = wahl_polarization(U236, boundary - 1e-8)
        above = wahl_polarization(U236, boundary + 1e-8)
        @test below ≈ above atol = 1e-6
    end

    # σ_Z is not, and that is the model rather than an implementation slip. Eqs. (12c) and (12g)
    # impose σ_Z(50) as a flat value over the two legs adjoining the Z = 50 crossing while the
    # crossing itself keeps the peak value, and eq. (14c) returns the far wing to σ_Z(B5) rather
    # than continuing the wing slope. Fig. 18a, p. 31, draws both as steps.
    @test wahl_dispersion(U236, U236.B4 - 1e-8) ≈ U236.σZ50 rtol = RTOL
    @test wahl_dispersion(U236, U236.B4 + 1e-8) > U236.σZ50
    @test wahl_dispersion(U236, U236.Bb - 1e-8) > U236.σZ50
    @test wahl_dispersion(U236, U236.B6 - 1e-8) < wahl_dispersion(U236, U236.B6 + 1e-8)
    @test wahl_dispersion(U236, U236.B6 + 1e-8) ≈ wahl_dispersion(U236, U236.B5) rtol = RTOL

    # The Bb–B4 leg is constructed so that it reaches ΔZ_max exactly at the crossing.
    @test wahl_polarization(U236, U236.Bb) ≈ U236.ΔZmax atol = 1e-6

    # In the peak region the distribution is the straight line of eq. (11a) about A' = 140.
    @test wahl_polarization(U236, 140) ≈ -0.487 rtol = RTOL
    @test wahl_polarization(U236, 150) ≈ -0.487 - 0.0080 * 10 rtol = RTOL
    @test wahl_dispersion(U236, 140) ≈ 0.566 rtol = RTOL
    @test wahl_dispersion(U236, 150) ≈ 0.566 - 0.0038 * 10 rtol = RTOL

    # σ_Z narrows to σ_Z(50) through the near-symmetry region, the signature of the shell
    # closure: the charge is pinned at Z = 50 and the distribution about it is tighter.
    @test wahl_dispersion(U236, (U236.Bb + U236.B4) / 2) ≈ U236.σZ50 rtol = RTOL
    @test U236.σZ50 < U236.σZ140

    # 252-Cf, a heavier and more charged system: every difference in eq. (17) contributes.
    CF252 = WahlSystematics(98, 252, 0.0)
    @test CF252.σZ140 ≈ 0.566 + 0.0064 * 16 + 0.0109 * (0 - 6.551) rtol = RTOL
    @test CF252.ΔZ140 ≈ -0.487 + 0.0180 * 16 - 0.00203 * 16^2 rtol = RTOL
    @test CF252.σZ50 ≈ 0.356 + 0.060 * 6 rtol = RTOL
    @test wahl_polarization(CF252, 126) ≈ 0.0 atol = 1e-12

    # Outside the fitted ranges the constructor refuses rather than extrapolating quietly.
    @test_throws ArgumentError WahlSystematics(88, 236, 6.5)      # Z_F below 90
    @test_throws ArgumentError WahlSystematics(92, 260, 6.5)      # A_F above 252
    @test_throws ArgumentError WahlSystematics(92, 236, 12.0)     # past the low-energy branch
    @test_throws ArgumentError WahlSystematics(92, 236, -1.0)
    # And the light fragment is not its domain: ΔZ there is the complement, not a lookup.
    @test_throws ArgumentError wahl_polarization(U236, 100)
    @test_throws ArgumentError wahl_dispersion(U236, 100)
end

@testset "where the systematics applies at all" begin
    table = mass_table()
    # Inside the fitted range: every shipped system is, which is what lets a run offer the
    # systematics beside whatever it used.
    @test is_wahl_applicable(table, spontaneous_fission(Nuclide(98, 252)))
    @test is_wahl_applicable(
        table,
        neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth"),
    )
    @test is_wahl_applicable(
        table,
        neutron_induced_fission(Nuclide(92, 233), 2.53e-8, "nth"),
    )
    @test is_wahl_applicable(
        table,
        neutron_induced_fission(Nuclide(94, 239), 2.53e-8, "nth"),
    )

    # Outside it, the predicate says so rather than leaving a caller to catch the constructor.
    # 232-Th is a trap here and not a case: Z_F = 90 and A_F = 232 are both inside the fitted
    # range, so it is applicable.
    @test is_wahl_applicable(table, spontaneous_fission(Nuclide(90, 232)))
    @test !is_wahl_applicable(table, spontaneous_fission(Nuclide(88, 226)))   # Z_F below 90
    @test !is_wahl_applicable(table, spontaneous_fission(Nuclide(98, 254)))   # A_F above 252
    @test !is_wahl_applicable(
        table,
        neutron_induced_fission(Nuclide(92, 235), 15.0, "nfast"),
    )

    # And it agrees with the constructor wherever the constructor has an opinion.
    for system in (
        spontaneous_fission(Nuclide(98, 252)),
        neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth"),
        neutron_induced_fission(Nuclide(92, 235), 15.0, "nfast"),
    )
        constructed = try
            WahlSystematics(table, system) !== nothing
        catch err
            err isa ArgumentError || rethrow()
            false
        end
        @test is_wahl_applicable(table, system) == constructed
    end
end

@testset "Wahl systematics as a charge model" begin
    table = mass_table()
    system = neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth")
    systematics = WahlSystematics(table, system)
    @test systematics isa ZpModel
    @test charge_polarization(systematics, 140) == wahl_polarization(systematics, 140)
    @test charge_dispersion(systematics, 140) == wahl_dispersion(systematics, 140)

    # The light side is the complement, eqs. (9c) and (9d).
    @test most_probable_charge(systematics, 100) ≈
          92 - most_probable_charge(systematics, 136) rtol = RTOL
    @test charge_dispersion(systematics, 100) == charge_dispersion(systematics, 136)

    # Reduced to a plain Gaussian per mass, the width is the lattice second moment, above σ_Z.
    effective = effective_charge_distribution(systematics, 118:160)
    @test effective isa ChargeDistribution
    @test charge_dispersion(effective, 140) > wahl_dispersion(systematics, 140)
    @test charge_polarization(effective, 300) == -0.5

    # The precursor excitation energy of a thermal capture is Sₙ of the compound nucleus, so the
    # systematics built from the system agrees with the one built from that number by hand.
    excitation = value(compound_nucleus_excitation(table, system))
    @test systematics.PE ≈ excitation rtol = RTOL
    @test WahlSystematics(92, 236, excitation).ΔZ140 ≈ systematics.ΔZ140 rtol = RTOL
    @test WahlSystematics(table, spontaneous_fission(Nuclide(98, 252))).PE == 0.0
    @test WahlSystematics(table, system; neutron_pairing = true).neutron_pairing
    @test_throws ArgumentError effective_charge_distribution(systematics, Int[])
end

guarded(CHARGE_DISTRIBUTION_AVAILABLE, "Wahl systematics against the evaluated tables") do
    @testset "Wahl systematics against the evaluated tables" begin
        # An external check, against numbers this package did not produce. The shipped tables are a
        # per-reaction least-squares fit of the same model — Fig. 18 quotes reduced χ² of 2.9 for
        # that against 7.9 for these systematics — so they agree in structure without coinciding, and
        # what is pinned here is structure and size, not identity. A transcription error in Table 2
        # would break all of it; the internal tests above would not notice.
        table = mass_table()
        systems = (
            (
                "235-U(n_th,f)",
                U235_CHARGE_DISTRIBUTION_FILE,
                neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth"),
                0.11,
            ),
            (
                "252-Cf(sf)",
                CF252_CHARGE_DISTRIBUTION_FILE,
                spontaneous_fission(Nuclide(98, 252)),
                0.30,
            ),
        )

        for (name, path, system, tolerance) in systems
            isfile(path) || continue
            evaluated = read_charge_distribution(path)
            systematics = WahlSystematics(table, system)

            # Mean absolute difference over the whole heavy branch, held loosely because the two
            # parameter sets are a systematics and a per-reaction fit of the same model and are
            # not meant to coincide.
            #
            # Where the difference sits is not the same in every case, and it is worth not
            # assuming. Measured over the shipped tables: for 252-Cf the Zₚ = 50 crossing
            # dominates — the largest difference is 0.75 at A = 130 against Bb = 126.8, and the
            # mean falls from 0.26 to 0.14 with the crossing region excluded — because SL50 is
            # extrapolated 16 mass units past the data it was fitted on and a small horizontal
            # offset costs a large pointwise difference on the steepest leg. For 235-U the
            # largest difference is also at the crossing but the mean barely moves without it,
            # 0.062 to 0.059: the disagreement is a broad offset across the branch. For 239-Pu
            # the crossing is not where it lives at all — the largest difference, 0.25, is at
            # A = 157, inside the peak-region leg, and excluding the crossing makes the mean
            # worse rather than better.
            differences =
                [abs(wahl_polarization(systematics, A) - ΔZ) for (A, ΔZ) in evaluated.ΔZ]
            @test sum(differences) / length(differences) < tolerance

            # The narrow σ_Z band is the crossing made visible. Both the model and the table must
            # have one, and they must agree on where it is to within a few mass units.
            masses = sort(collect(keys(evaluated.σ_Z)))
            widths = [evaluated.σ_Z[A] for A in masses]
            threshold = minimum(widths) + 0.25 * (maximum(widths) - minimum(widths))
            tabulated_band = [A for (A, w) in zip(masses, widths) if w < threshold]
            @test !isempty(tabulated_band)
            @test abs(minimum(tabulated_band) - systematics.Bb) < 3.0
            @test abs(
                (minimum(tabulated_band) + maximum(tabulated_band)) / 2 -
                (systematics.Bb + systematics.B4) / 2,
            ) < 10.0

            # Bb is the mass at which ΔZ attains ΔZ_max, so the table's own maximum is where it
            # should be. This is what settles eq. (12): the blended A'_max lands within a mass unit
            # of it, and either estimate taken alone does worse. Zₚ there is recorded alongside,
            # since it is the quantity the interval is named for and the blend does not reach 50
            # exactly except where F₁ vanishes.
            polarizations = [evaluated.ΔZ[A] for A in masses]
            @test abs(masses[argmax(polarizations)] - systematics.Bb) < 1.5
            crossing =
                systematics.Bb * systematics.Z_F / systematics.A_F + systematics.ΔZmax
            @test 49.0 < crossing <= 50.0 + 1e-9
        end

        # Where A' = 140 falls inside the peak region the model is at its best determined, and the
        # agreement is then a real number rather than a bound. 240-Pu is the case: σ_Z(140) is
        # 0.59141 against a tabulated 0.59131.
        plutonium = joinpath(DATA, "Pu239_nth", "charge_distribution_vs_A.dat")
        if isfile(plutonium)
            evaluated = read_charge_distribution(plutonium)
            systematics = WahlSystematics(
                table,
                neutron_induced_fission(Nuclide(94, 239), 2.53e-8, "nth"),
            )
            @test 140 > systematics.B4
            @test wahl_dispersion(systematics, 140) ≈ charge_dispersion(evaluated, 140) atol =
                1e-3
        end
    end
end

@testset "even-odd factors" begin
    table = mass_table()
    system = neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth")
    w = WahlSystematics(table, system)

    # Eq. (12a): F(A) = 1 through the near-symmetry region, so both factors are unity there.
    @test wahl_even_odd_factors(w, (w.Ba + w.Bb) / 2) == (1.0, 1.0)
    @test wahl_even_odd_factors(w, w.B4 - 1e-8) == (1.0, 1.0)

    # Eqs. (11d), (11e): constant through the peak region, at the Table 2 values.
    peak = wahl_even_odd_factors(w, 140)
    @test peak == wahl_even_odd_factors(w, w.B4 + 1e-8)
    @test peak[1] ≈ 1.207 rtol = RTOL
    @test peak[2] ≈ 1.076 rtol = RTOL
    # Both exceed unity, which is what an even-odd enhancement means.
    @test all(>(1), peak)

    # Eqs. (13f), (13h): the wing carries a slope, the far wing returns to the peak value.
    @test wahl_even_odd_factors(w, w.B5 + 5)[1] > peak[1]
    @test wahl_even_odd_factors(w, w.B6 + 1e-8) == peak

    # Eq. (9): the Gaussian integrated over unit charge intervals, modulated by F and
    # renormalized, so the yields sum to one, and the modulation follows the parity of Z. The
    # product-level form carries both factors; at A even the parities of Z and N move together.
    numbers, yields = fractional_independent_yields(w, 140, 0.0)
    @test sum(yields) ≈ 1.0 atol = 1e-12
    @test all(>=(0), yields)
    Zₚ = most_probable_charge(w, 140)
    σ = wahl_dispersion(w, 140)
    plain = [erf_interval(Z, Zₚ, σ) for Z in numbers]
    plain ./= sum(plain)
    ratio = yields ./ plain
    inside = [i for i in eachindex(numbers) if plain[i] > 1e-6]
    @test all(ratio[i] > 1 for i in inside if iseven(numbers[i]))
    @test all(ratio[i] < 1 for i in inside if isodd(numbers[i]))
    # And the lattice form, not the density: the second moment is σ² + 1/12 without the factors.
    flat = [erf_interval(Z, 54.3, 0.566) for Z in 40:70]
    flat ./= sum(flat)
    @test charge_moments(collect(40:70), flat)[2] ≈ sqrt(0.566^2 + 1 / 12) rtol = 1e-3
    @test_throws ArgumentError wahl_even_odd_factors(w, 100)
end

guarded(
    CHARGE_DISTRIBUTION_AVAILABLE && U235_POLARIZATION_AVAILABLE,
    "the 1988 model is the closest of the layers to the supplied tables",
) do
    @testset "the 1988 model is the closest of the layers to the supplied tables" begin
        # The supplied tables are of the same lineage as Wahl (1988) but not a reproduction of
        # it. Reduced to effective Gaussians, the 1988 model at fragment level lies within a few
        # hundredths of them in both ΔZ and width, and nearer than the systematics everywhere.
        table = mass_table()
        for (path, system, bound) in (
            (
                U235_CHARGE_DISTRIBUTION_FILE,
                neutron_induced_fission(Nuclide(92, 235), 2.53e-8, "nth"),
                0.035,
            ),
            (CF252_CHARGE_DISTRIBUTION_FILE, spontaneous_fission(Nuclide(98, 252)), 0.025),
        )
            evaluated = read_charge_distribution(path)
            masses = sort([A for A in keys(evaluated.σ_Z) if 2A >= system.compound.A])
            reaction = effective_charge_distribution(charge_model(table, system), masses)
            systematics =
                effective_charge_distribution(WahlSystematics(table, system), masses)
            mad(d, f) = sum(abs(f(d, A) - f(evaluated, A)) for A in masses) / length(masses)
            @test mad(reaction, charge_polarization) < bound
            @test mad(reaction, charge_dispersion) < bound
            @test mad(reaction, charge_polarization) < mad(systematics, charge_polarization)
            @test mad(reaction, charge_dispersion) < mad(systematics, charge_dispersion)
        end
    end
end
