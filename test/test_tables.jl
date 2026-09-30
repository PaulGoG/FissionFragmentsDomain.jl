@testset "table specs" begin
    spec = TableSpec([:A, :yield, :uncertainty]; skip = 1, optional = (:uncertainty,))
    @test spec.columns == [:A, :yield, :uncertainty]
    @test spec.skip == 1

    @test_throws ArgumentError TableSpec(Symbol[])
    @test_throws ArgumentError TableSpec([:A, :A])
    @test_throws ArgumentError TableSpec([:A, :B]; optional = (:C,))
    # Optional columns have to form a suffix, or "a leading prefix of the optional ones" is
    # not a well-defined file shape.
    @test_throws ArgumentError TableSpec([:A, :B, :C]; optional = (:A,))
    @test_throws ArgumentError TableSpec([:A]; skip = -1)
end

@testset "reading tables" begin
    directory = mktempdir()

    # The shape ExforFissionData.jl emits when the measurements carry uncertainties.
    with_uncertainty = joinpath(directory, "with.dat")
    write(with_uncertainty, "A yield erryield\n100 6.0 0.1\n101 6.5 0.2\n")

    # The shape it emits when they do not — the same reader must accept both, which is what
    # retires the prototype's separate column-padding script.
    without = joinpath(directory, "without.dat")
    write(without, "A yield\n100 6.0\n101 6.5\n")

    spec = TableSpec([:A, :yield, :uncertainty]; skip = 1, optional = (:uncertainty,))

    table = read_delimited_table(with_uncertainty, spec)
    @test length(table) == 2
    @test has_column(table, :uncertainty)
    @test column(table, :A, Int) == [100, 101]
    @test column(table, :yield, Float64) == [6.0, 6.5]
    @test column_or(table, :uncertainty, Float64, 0.0) == [0.1, 0.2]

    bare = read_delimited_table(without, spec)
    @test !has_column(bare, :uncertainty)
    @test column_or(bare, :uncertainty, Float64, 0.0) == [0.0, 0.0]
    @test_throws ArgumentError column(bare, :uncertainty, Float64)

    # Repeated delimiters collapse, so the space-padded legacy layouts need no separate path.
    padded = joinpath(directory, "padded.dat")
    write(padded, "Z A symbol excess uncertainty\n   0    1    n     8071.3    0.00044\n")
    padded_table = read_delimited_table(
        padded,
        TableSpec([:Z, :A, :symbol, :excess, :uncertainty]; skip = 1),
    )
    @test column(padded_table, :excess, Float64) == [8071.3]

    # A ragged row is named by file and line rather than parsed into silence.
    ragged = joinpath(directory, "ragged.dat")
    write(ragged, "A yield\n100 6.0\n101\n")
    @test_throws ArgumentError read_delimited_table(
        ragged,
        TableSpec([:A, :yield]; skip = 1),
    )

    empty_file = joinpath(directory, "empty.dat")
    write(empty_file, "A yield\n")
    @test_throws ArgumentError read_delimited_table(
        empty_file,
        TableSpec([:A, :yield]; skip = 1),
    )

    @test_throws ArgumentError read_delimited_table(joinpath(directory, "absent.dat"), spec)

    unparseable = joinpath(directory, "bad.dat")
    write(unparseable, "A yield\n100 six\n")
    bad = read_delimited_table(unparseable, TableSpec([:A, :yield]; skip = 1))
    exception = try
        column(bad, :yield, Float64)
        nothing
    catch error
        error
    end
    @test exception isa ArgumentError
    @test occursin("row 1", exception.msg)
    @test occursin("yield", exception.msg)
end
