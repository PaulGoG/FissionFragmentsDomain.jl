@testset "jointed segments" begin
    segments = SegmentedCurve([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)])
    @test length(segments) == 4

    # Endpoints are reproduced exactly, and interpolation is linear between them.
    @test segments(126) ≈ 0.5 rtol = RTOL
    @test segments(130) ≈ 0.14 rtol = RTOL
    @test segments(150) ≈ 0.607 rtol = RTOL
    @test segments(174) ≈ 0.865 rtol = RTOL
    @test segments(140) ≈ 0.14 + (0.607 - 0.14) * (140 - 130) / (150 - 130) rtol = RTOL

    # Undefined above the last point; the published parameterizations end at A_H_max.
    @test segments(175) === nothing

    # Undefined below the first point as well: a fit is not extrapolated on either side.
    @test segments(125) === nothing

    # A segment crossing zero leaves the function undefined there rather than negative.
    crossing = SegmentedCurve([(120, 1.0), (140, -1.0)])
    @test crossing(130) === nothing
    @test crossing(121) ≈ 0.9 rtol = RTOL

    @test_throws ArgumentError SegmentedCurve([(126, 0.5)])
    @test_throws ArgumentError SegmentedCurve([(150, 0.5), (126, 0.1)])
    @test_throws ArgumentError SegmentedCurve([(126, 0.5), (126, 0.1)])
end

@testset "reading a tabulated curve" begin
    # The layout a temperature-ratio run writes: a header line, then A_H, R_T and an uncertainty
    # column that the partition discards because it carries none downstream.
    mktempdir() do dir
        path = joinpath(dir, "R_T_vs_A_H_segmented_systematic_trend_run.csv")
        write(
            path,
            """
            A_H,R_T,R_T_uncertainty
            126,1.0,0.0
            127,1.205662,0.007435
            128,1.471892,0.02157
            129,1.932797,0.061908
            """,
        )
        segments = read_segmented_curve(path)
        @test length(segments) == 4
        @test segments.source == path

        # Tabulated at every integer mass, so evaluation reproduces each value exactly rather
        # than interpolating across it.
        @test segments(126) == 1.0
        @test segments(127) == 1.205662
        @test segments(128) == 1.471892
        @test segments(129) == 1.932797

        # The range is the fit's, not the configured domain's: a set that stops short leaves
        # the partition undefined above its last point instead of extrapolating.
        @test segments(130) === nothing

        # The uncertainty column is optional, as it is for every other reader here.
        bare = joinpath(dir, "bare.csv")
        write(bare, "A_H,R_T\n126,1.0\n127,1.2\n")
        @test read_segmented_curve(bare)(127) == 1.2

        # Columns are taken by position, so a header naming something else entirely reads the
        # same. Nothing downstream may be coupled to the text of a header.
        renamed = joinpath(dir, "renamed.csv")
        write(renamed, "heavy_mass,ratio,sigma\n126,1.0,0.0\n127,1.2,0.0\n")
        @test read_segmented_curve(renamed)(127) == 1.2

        # The excitation-ratio spec reads the same shape under its own quantity name.
        excitation = joinpath(dir, "excitation.csv")
        write(excitation, "A_H,excitation_ratio\n126,0.5\n174,0.865\n")
        @test read_segmented_curve(excitation; spec = EXCITATION_RATIO_CURVE_SPEC)(126) ==
              0.5
    end
end
