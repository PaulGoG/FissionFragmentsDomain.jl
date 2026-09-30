# Data

Two tables ship with this package: `reference/mass_excess_ame2020.dat`, reached as
`AME2020_MASS_EXCESS_FILE`, and `reference/shell_corrections_gilbert_cameron_1965.dat`, reached as
`GILBERT_CAMERON_SHELL_CORRECTION_FILE`. Every other table the readers accept is third-party data
that this package does not redistribute. It belongs in the data directory of the calculation that consumes
it, laid out as below; the readers take any path.

## Layout

```
data/
├── reference/                      evaluations that belong to no one system
│   ├── mass_excess_ame2020.dat     shipped
│   └── shell_corrections_gilbert_cameron_1965.dat   shipped
└── <system>/                       Cf252_sf, U235_nth, …, the token of system_label
    ├── charge_distribution_vs_A.dat
    ├── Y_vs_A_TKE.dat
    ├── Y_vs_A/*.dat, TKE_vs_A/*.dat
    └── <temperature-ratio run>/    manifest_<run>.toml and the files it names
```

| Reader | Layout, one header row | Units |
|---|---|---|
| `read_mass_excess_table` | `Z A element mass_excess mass_excess_uncertainty` | keV |
| `read_shell_correction_table` | `n S_N S_Z` | MeV |
| `read_charge_distribution` | `A dZ sigma_Z` | charge units |
| `read_mass_energy_yield` | `A TKE Y [Y_uncertainty]` | MeV, % |
| `read_mass_yield`, `factorized_yield` | `A Y [Y_uncertainty]` | % |
| `factorized_yield` | `A TKE [TKE_uncertainty]` | MeV |
| `read_kinetic_energy_dispersion` | `A sigma_TKE` | MeV |
| `read_segmented_curve` | `A_H,R_T[,R_T_uncertainty]` or `A_H,excitation_ratio[,…]`, comma-separated | — |
| `read_temperature_ratio_manifest` | TOML, below | — |

A file name is `<quantity>_vs_<abscissa>`. Where a quantity was measured more than once, the
directory carries the quantity and the file carries the provenance,
`<accession>_<Author>_<year>.dat`, as `ExforFissionData.jl` writes it.

All files are whitespace-separated with a single header line, except the `.csv` tables, which are
comma-separated. **Readers take columns by position, not by header text**, and the uncertainty
column may be omitted rather than zero-filled. Renaming a header is therefore a no-op for code,
which is what made this layout safe to adopt.

Nucleon-number columns are read through `integer_column`, which accepts an integral value
written in floating point — some of these tables were emitted by a generator and write `118.0`
where they mean `118`.

## The shipped mass table

**`reference/mass_excess_ame2020.dat`** — the 2020 atomic mass evaluation.
W. J. Huang, M. Wang, F. G. Kondev, G. Audi, S. Naimi, *Chinese Physics C* **45**, 030002 (2021),
[doi:10.1088/1674-1137/abddb0](https://doi.org/10.1088/1674-1137/abddb0);
M. Wang, W. J. Huang, F. G. Kondev, G. Audi, S. Naimi, *Chinese Physics C* **45**, 030003 (2021),
[doi:10.1088/1674-1137/abddaf](https://doi.org/10.1088/1674-1137/abddaf).
The source is the unrounded mass table `mass.mas20` of the Atomic Mass Data Center, published at
<https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt> (file dated 3 March 2021, MD5
`6d28b75833cf53c7cc230223f63da6f6`). The shipped table carries that file's mass excess and
uncertainty for all 3558 nuclides, including the 1008 values the evaluation marks with `#` as
estimated — the mark itself is not carried, so the table does not tell an estimated mass from a
measured one — after one pass through single precision: every entry is the source rounded to
the nearest `Float32` and printed with eight decimals, verified row by row against the source on
23 September 2026. The pass moves a mass excess by at most 0.006 keV, and by at most 0.004 keV
over 70 ≤ A ≤ 180, where it exceeds the evaluation's own uncertainty for ⁸⁶Kr alone, by 0.04 eV;
elsewhere only for ¹H, ²H, ³H and ³He. No quantity computed here is sensitive to it, and the
table is kept as it is so that the reference values the tests pin keep their meaning.
`HEADERLESS_MASS_EXCESS_SPEC` reads the same five columns without a header line, the form the
table circulated in before the shipped copy was given one.

## Shell corrections

**`reference/shell_corrections_gilbert_cameron.dat`** — shell corrections `S(N)` and `S(Z)`,
Table III of A. Gilbert and A. G. W. Cameron, *Canadian Journal of Physics* **43**, 1446 (1965),
pp. 1453–1455, as redistributed in the IAEA Reference Input Parameter Library. The paper
tabulates `S(Z)` only to `Z = 98`; the file is rectangular and pads the column with zeros beyond
that, which the reader discards rather than storing as a correction of zero.

## Charge distributions

No charge-distribution table ships, and none is needed. `charge_model` builds the isobaric charge
distribution from Wahl's `Zₚ` model, whose every number traces to a public source:

- **The reaction's own parameters.** Table A of A. C. Wahl, *At. Data Nucl. Data Tables* **39**, 1
  (1988), doi:10.1016/0092-640X(88)90016-2, fits them for ²³⁵U(n_th,f), ²³³U(n_th,f),
  ²³⁹Pu(n_th,f) and ²⁵²Cf(sf). `Wahl1988` reproduces that evaluation's own Tables I–IV row by
  row, and the test suite checks it against them, transcribed in
  `test/references/wahl1988/`.
- **The systematics.** Eq. (17) and Table 2 of A. C. Wahl, *Systematics of Fission-Product Yields*,
  LA-13928 (2002), doi:10.2172/809574, estimate the parameters from `Z_F`, `A_F` and the
  excitation energy, for any reaction in `90 ≤ Z_F ≤ 98`, `230 ≤ A_F ≤ 252`, `PE ≤ 8 MeV`.
- **The conventional means**, `ΔZ = −0.5` and `σ_Z = 0.6`, elsewhere.

Fig. 18 of LA-13928 quotes a reduced `χ²` of 2.9 for the first against 7.9 for the second, on
²³⁵U(n_th,f).

A consumer holding an evaluated table of its own can still read it with
`read_charge_distribution`, in the `A dZ sigma_Z` layout above. Such a table is an input of the
consumer's, not of this package.

## Temperature-ratio run records

A temperature-ratio extraction writes the whole results directory of one run, and a consumer
stages all of it, not one table out of it:

```
results/<system>/
├── manifest_<run>.toml                          the run record, read first
├── R_T_vs_A_H_segmented_<label>_<run>.csv       header: A_H,R_T,R_T_uncertainty
├── r_nu_vs_A_H_pivots_<label>_<run>.csv         header: A_H,r_nu,r_nu_uncertainty
├── total_average_R_T_<run>.csv
├── dataset_diagnostics_<run>.csv
└── metadata_<run>.toml
```

`staged_manifest` finds the record in the directory it was copied to, and a consumer picks one of
its `[[segmented_curve]]` entries by `label`; the systematic trend carries
`SYSTEMATIC_TREND_LABEL`. Keep one staged run per directory — two manifests in one directory are
refused, because nothing would say which run to take.

Copying one table by hand is what this replaces, and the reason is physical rather than tidy: the
two CSVs above hold different quantities in the same column position, `R_T ≈ 1.2` against
`r_ν ≈ 0.5`, and both are plausible. Staging the pivots file by mistake runs to completion and
reports a number that looks like a result. The manifest is what makes the quantity checkable —
`read_temperature_ratio_manifest` refuses a manifest that does not declare `ordinate = "R_T"` and
`abscissa = ["A_H"]`, and resolves every curve through `temperature_ratio_file` and through no
other field. The `columns` the manifest states are documentation: tables are read by position.

The `[system]` table is `system_record` of the system the extraction ran for. The level density
parameter ratio that produced `R_T` must be averaged over the charge distribution in the same
order as the partition that applies it; see `RatioAveraging`.

## Gilbert–Cameron shell corrections

`reference/shell_corrections_gilbert_cameron_1965.dat` is Table III, "Pairing energies and shell
corrections", of A. Gilbert and A. G. W. Cameron, *Can. J. Phys.* **43**, 1446 (1965),
pp. 1453–1455, doi:10.1139/p65-139. It holds the shell corrections `S(N)` for `N = 11–150` and
`S(Z)` for `Z = 11–98`, in MeV, with the paper's two decimals. The paper tabulates no `S(Z)` past
`Z = 98`, and those cells read `NaN`. The pairing energies `P(Z)`, `P(N)` of the same table are
not used by the level density parameter and are not shipped.

The corrections are the residuals of the exponential semi-empirical mass formula of A. G. W.
Cameron and R. M. Elkin, *Can. J. Phys.* **43**, 1288 (1965), doi:10.1139/p65-123, split into a
function of `Z` and a function of `N`. The paper states them good to about 200 keV
(p. 1455).

The table was transcribed three times independently from page images of the printed article,
reading row by row, column by column, and even rows before odd ones. The three agree in every
cell. They also agree in all 228 values with an independent electronic copy, the tables of
Geant4's `G4CameronGilbertShellCorrections.cc` (github.com/Geant4/geant4, commit `f3d5293d38`). Checks: `Σ S(N) = 1496.25 MeV` over 140 values, `Σ S(Z) = −1018.24 MeV` over 88 values, and
the closed shells `Z = 28, 50, 82` and `N = 28, 50, 82, 126` are local minima. A copy that pads
`S(Z)` past `Z = 98` with `0.00` reads to the same table.
