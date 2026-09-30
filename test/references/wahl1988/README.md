# Wahl (1988), Tables I–IV

Columns of Tables I–IV of A. C. Wahl, "Nuclear-charge distribution and delayed-neutron yields
for thermal-neutron-induced fission of ²³⁵U, ²³³U, and ²³⁹Pu and for spontaneous fission of
²⁵²Cf", *At. Data Nucl. Data Tables* **39**, 1–156 (1988), doi:10.1016/0092-640X(88)90016-2,
pp. 54–61. There is one file per reaction, named by Wahl's reaction symbol, with one row per
product mass number `A`.

| Column | Wahl's heading | Meaning |
|---|---|---|
| `A` | A | product mass number, after prompt-neutron emission |
| `NU` | NU(A) | `ν̄_A`, prompt neutrons emitted to form products of mass `A`; `A' = A + ν̄_A` |
| `Zp` | Zp(A) | most probable charge of the `Zₚ` model |
| `AVEZ_calc`, `RMS_calc` | AVE.Z(A), RMS(A) CALC | mean charge and rms dispersion of the `Zₚ`-model yields |
| `AVEZ_expt`, `RMS_expt` | AVE.Z(A), RMS(A) EXPT | the same from the evaluated experimental yields; `NaN` where Wahl gives none |

The values were transcribed from the scanned paper three times, independently. Two copies were
made from page renders at about 100 dpi and one from 300-dpi strips. Every cell on which the
three disagreed, 15 of about 4800, was read again from the 300-dpi scan. The values are as
printed; nothing is corrected.

These are reference values for the tests, not inputs. `test_wahl1988.jl` requires the model to
reproduce the calculated columns from `A` and `NU` alone.
