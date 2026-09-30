# Reference data

Output of the implementation submitted with the 2023 MSc thesis, kept here so this package can
be pinned against it. The configuration these files come from:

| | |
|---|---|
| Fissioning system | ²⁵²Cf(SF) |
| Heavy fragment masses | 126 … 174 |
| TKE grid | 130 … 230 MeV, step 2 MeV |
| Charge numbers per mass | 5 |
| Mass excesses | the 2020 atomic mass evaluation |
| Charge distribution | the per-reaction `ΔZ(A)`, `σ_Z(A)` fit for ²⁵²Cf |
| Evaporation cross section | variable |
| Level density parameter | back-shifted Fermi gas |

The quantities kept here precede the TXE partition, so they are the same in both archived runs.
Files are named by the `<quantity>_vs_<abscissa>` rule used for every table in the toolchain.

| File, under `Cf252_sf/` | Content |
|---|---|
| `Q_vs_A_H_Z_H.dat` | `Q(A_H, Z_H)` in MeV, 245 rows — one per fragmentation |
| `a_vs_A_H_Z_H.dat` | `a(A, Z)` in MeV⁻¹ over the fragment domain |
| `Y_vs_A.dat`, `Y_vs_Z.dat`, `Y_vs_TKE.dat`, `TKE_vs_A_H.dat` | the fragment yield projections |

Only the header lines were rewritten when the naming convention was applied; every data row is
byte-identical to what the archived run produced, and `.gitattributes` keeps these files exempt
from line-ending normalisation so they stay that way.

These are not inputs. They exist so that a change in the physics shows up as a failing test
rather than as a number nobody checked.
