"""
    FissionFragmentsDomain

The fragments of binary fission, as the models of prompt emission and of the excitation-energy
partition see them.

A fissioning system splits into complementary pairs `(A_L, Z_L) + (A_H, Z_H)` at a total kinetic
energy `TKE`. This package holds what those models share: nuclides and fissioning systems, the
atomic mass evaluation and the `Q`-values it gives, the isobaric charge distribution (evaluated
tables, or Wahl's `Zₚ` model with a reaction's own 1988 parameters or the 2002 systematics),
the fragmentation domain built from it, fragment level
density parameters and their ratio across a pair, the energy balance `TXE = Q + E*_CN − TKE`, the
relation between the temperature ratio `R_T` and the excitation-energy partition, the fragment
yield `Y(A, Z, TKE)` and averages over it, and the heavy-fragment curves `R_T(A_H)` and
`E*_H/TXE (A_H)` together with the run record that carries a temperature-ratio extraction to its
consumers.
"""
module FissionFragmentsDomain

using Measurements: Measurement, measurement, value
using RelocatableFolders: @path
using SpecialFunctions: erf
using TOML: TOML

export Nuclide, NEUTRON, neutron_number, element_symbol
export FissioningSystem,
    CHANNELS,
    CHANNEL_REACTION,
    REACTIONS,
    spontaneous_fission,
    neutron_induced_fission,
    is_spontaneous,
    reaction,
    symmetric_mass,
    system_label,
    system_notation,
    system_record
export TableSpec,
    DelimitedTable, read_delimited_table, column, column_or, integer_column, has_column
export AME2020_MASS_EXCESS_FILE,
    MassExcessTable,
    MASS_EXCESS_SPEC,
    HEADERLESS_MASS_EXCESS_SPEC,
    read_mass_excess_table,
    mass_excess,
    separation_energy,
    neutron_separation_energy,
    q_value,
    complementary_fragment
export ChargeModel,
    ChargeDistribution,
    CHARGE_DISTRIBUTION_SPEC,
    DEFAULT_CHARGE_POLARIZATION,
    DEFAULT_CHARGE_DISPERSION,
    read_charge_distribution,
    mean_charge_distribution,
    charge_polarization,
    charge_dispersion,
    unchanged_charge_distribution,
    most_probable_charge,
    charge_probability,
    charge_numbers,
    fragment_most_probable_charge,
    fragment_charge_probability,
    charge_model_label
export ZpModel,
    fractional_independent_yields,
    fragment_charge_yields,
    charge_moments,
    even_odd_factors,
    effective_charge_distribution,
    Wahl1988,
    WAHL_1988,
    WAHL_1988_TABLE_A,
    wahl_1988,
    charge_model
export WahlSystematics,
    WAHL_LOW_ENERGY,
    WAHL_VALIDITY,
    is_wahl_applicable,
    wahl_polarization,
    wahl_dispersion,
    wahl_even_odd_factors
export Fragmentation,
    FragmentationDomain,
    fragmentation_domain,
    fragments,
    distinct_splits,
    mass_range,
    is_symmetric,
    is_self_complementary,
    split_weight,
    charges,
    symmetric_charge_set_is_invariant
export Energetics,
    Exclusion,
    energetics,
    sweep_energetics,
    exclusion_summary,
    compound_nucleus_excitation,
    total_excitation_energy,
    kinetic_energy,
    fragment_energy_per_nucleon
export LevelDensityModel,
    LEVEL_DENSITY_MODELS,
    GILBERT_CAMERON_SHELL_CORRECTION_FILE,
    ShellCorrectionTable,
    read_shell_correction_table,
    GilbertCameron,
    is_deformed_gilbert_cameron,
    BackShiftedFermiGas,
    level_density_parameter,
    shell_correction,
    liquid_drop_energy
export RatioAveraging,
    ChargeResolved,
    mean_total_excitation,
    RatioOfMeans,
    MeanOfRatios,
    level_density_ratio,
    heavy_excitation_fraction,
    temperature_ratio,
    RATIO_AVERAGINGS,
    ratio_averaging,
    ratio_averaging_label,
    level_density_label,
    temperature_ratio_slope
export KineticEnergyDispersion,
    KINETIC_ENERGY_DISPERSION_SPEC,
    DEFAULT_KINETIC_ENERGY_DISPERSION,
    read_kinetic_energy_dispersion,
    uniform_kinetic_energy_dispersion,
    kinetic_energy_dispersion,
    MassEnergyYield,
    MASS_ENERGY_YIELD_SPEC,
    MASS_YIELD_SPEC,
    MEAN_KINETIC_ENERGY_SPEC,
    read_mass_energy_yield,
    factorized_yield,
    symmetrized_yield,
    RECOMMENDED_MEAN_TOTAL_KINETIC_ENERGY,
    recommended_mean_total_kinetic_energy,
    fragment_charge_probabilities,
    FragmentYield,
    DEFAULT_YIELD_TOTAL,
    build_fragment_yield,
    coverage,
    yield_over,
    mean_kinetic_energy_by_mass,
    MassYield,
    mass_yield,
    read_mass_yield
export FragmentQuantity, average, total_average, weights_from
export SegmentedCurve,
    TEMPERATURE_RATIO_CURVE_SPEC, EXCITATION_RATIO_CURVE_SPEC, read_segmented_curve
export ManifestCurve,
    ManifestSystem,
    ManifestDomain,
    TemperatureRatioManifest,
    MANIFEST_CURVE_KINDS,
    MANIFEST_ORDINATE,
    MANIFEST_ABSCISSA,
    SYSTEMATIC_TREND_LABEL,
    read_temperature_ratio_manifest,
    write_temperature_ratio_manifest,
    staged_manifest,
    curve_labels,
    manifest_curve,
    temperature_ratio_path

"""
    AME2020_MASS_EXCESS_FILE

The shipped mass-excess table: the 2020 atomic mass evaluation, 3558 nuclides, in keV. Its
provenance is in `data/README.md`. The path survives relocation of the package, so it holds in
a system image as well as in a depot.

# Examples

```jldoctest
julia> length(read_mass_excess_table(AME2020_MASS_EXCESS_FILE))
3558
```
"""
const AME2020_MASS_EXCESS_FILE =
    @path joinpath(@__DIR__, "..", "data", "reference", "mass_excess_ame2020.dat")

"""
    GILBERT_CAMERON_SHELL_CORRECTION_FILE

The shipped shell corrections of Gilbert and Cameron, Table III of *Can. J. Phys.* **43**, 1446
(1965), pp. 1453–1455, doi:10.1139/p65-139: `S(N)` for `N = 11–150` and `S(Z)` for
`Z = 11–98`, in MeV, as printed. Rows past `Z = 98`, where the paper tabulates no `S(Z)`, carry
`NaN` in that column. The transcription and its checks, among them agreement with Geant4's
electronic copy in all 228 values, are in `data/README.md`. The path survives
relocation of the package, like [`AME2020_MASS_EXCESS_FILE`](@ref).

# Examples

```jldoctest
julia> corrections = read_shell_correction_table(GILBERT_CAMERON_SHELL_CORRECTION_FILE);

julia> corrections.S_Z[50], corrections.S_N[82]
(-19.83, 9.09)
```
"""
const GILBERT_CAMERON_SHELL_CORRECTION_FILE = @path joinpath(
    @__DIR__,
    "..",
    "data",
    "reference",
    "shell_corrections_gilbert_cameron_1965.dat",
)

# The version of this package, read once when it is compiled, for the run records it writes.
const PACKAGE_VERSION =
    VersionNumber(TOML.parsefile(joinpath(dirname(@__DIR__), "Project.toml"))["version"])

include("nuclides.jl")
include("tables.jl")
include("masses.jl")
include("charge.jl")
include("zp.jl")
include("domain.jl")
include("leveldensity.jl")
include("energetics.jl")
include("temperature_ratio.jl")
include("wahl.jl")
include("wahl1988.jl")
include("yields.jl")
include("averaging.jl")
include("curves.jl")
include("manifest.jl")

end
