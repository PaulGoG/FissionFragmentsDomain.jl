# The configuration of the archived reference run, and helpers shared by the test files.

const CF252 = spontaneous_fission(Nuclide(98, 252))
const CF252_HEAVY_MASSES = 126:174
const CF252_TKE = 130.0:2.0:230.0
const CF252_CHARGES_PER_MASS = 5

# Every fixture rests on a public source: the shipped AME2020 table, and for the charge
# distribution Wahl's 1988 model of ²⁵²Cf(sf), which reproduces that evaluation's own tables.
# The Gilbert–Cameron shell corrections are not shipped; the tests that need them read a data
# directory laid out as `data/README.md` describes, taken from `FISSION_FRAGMENTS_DOMAIN_TEST_DATA`
# when set, and are skipped when the file is absent.
const MASS_EXCESS_FILE = String(AME2020_MASS_EXCESS_FILE)
const SHELL_CORRECTION_FILE =
    joinpath(DATA, "reference", "shell_corrections_gilbert_cameron.dat")
const SHELL_CORRECTIONS_AVAILABLE = isfile(SHELL_CORRECTION_FILE)

# Relative tolerance of `isapprox` between Float64 operands when none is given.
const RTOL = sqrt(eps(Float64))

# Names of the testsets not run because their input data is absent.
const SKIPPED = String[]

"""
    guarded(flag, name, body)

Run `body()` when `flag` holds; otherwise record `name` in `SKIPPED`, so that a testset gated on
absent input data is reported at the end of the run rather than omitted without trace.

The method taking `body` first serves the `do`-block form, `guarded(flag, name) do … end`.
"""
function guarded(flag::Bool, name::AbstractString, body::Function)
    flag ? body() : push!(SKIPPED, String(name))
    return nothing
end

guarded(body::Function, flag::Bool, name::AbstractString) = guarded(flag, name, body)

mass_table() = read_mass_excess_table(MASS_EXCESS_FILE)

cf252_charge_distribution() = WAHL_1988[(98, 252)]

cf252_domain() = fragmentation_domain(
    CF252,
    cf252_charge_distribution(),
    CF252_HEAVY_MASSES;
    charges_per_mass = CF252_CHARGES_PER_MASS,
)

"""
    fragmentation_at(A_H, Z_H) -> Fragmentation

The ²⁵²Cf fragmentation with this heavy fragment, built without the charge distribution.

The reference tests pin `Q`, `a`, `E*_L`, `E*_H` and the emission sequences, none of which
depends on the isobaric charge probability — it is a yield weight, and these are per
configuration. Taking the fragmentation straight from `(A_H, Z_H)` keeps those tests running
without the unshipped systematics, and isolates them from the domain construction, which is
tested separately.
"""
function fragmentation_at(A_H::Integer, Z_H::Integer)
    heavy = Nuclide(Z_H, A_H)
    return Fragmentation(heavy, complementary_fragment(CF252, heavy), 1.0)
end

"""
    read_reference(name, columns) -> Dict

Read a reference file into a dictionary keyed by its leading integer columns.
"""
function read_reference(name::AbstractString, columns::Vector{Symbol}, keys::Int)
    table = read_delimited_table(joinpath(REFERENCE, name), TableSpec(columns; skip = 1))
    key_columns = [column(table, columns[i], Int) for i in 1:keys]
    values = column(table, columns[end], Float64)
    reference = Dict{NTuple{keys, Int}, Float64}()
    for row in eachindex(values)
        reference[ntuple(i -> key_columns[i][row], keys)] = values[row]
    end
    return reference
end

"The unit-interval integral of a normalised Gaussian, the lattice form of Wahl's eq. (7)."
erf_interval(Z, Zₚ, σ) =
    0.5 * (erf((Z - Zₚ + 0.5) / (σ * sqrt(2))) - erf((Z - Zₚ - 0.5) / (σ * sqrt(2))))
