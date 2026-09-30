# Reference data

## `Cf252_sf/`

Two tables over the ²⁵²Cf(sf) fragment domain, `A_H` 126–174, five charge numbers per mass. The
2023 MSc thesis implementation computed them first, and this package is pinned against them:

| File | Content | Rests on |
|---|---|---|
| `Q_vs_A_H_Z_H.dat` | `Q(A_H, Z_H)` in MeV, 245 rows, one per fragmentation | the 2020 atomic mass evaluation, doi:10.1088/1674-1137/abddb0 and doi:10.1088/1674-1137/abddaf |
| `a_vs_A_H_Z_H.dat` | `a(A, Z)` in MeV⁻¹ over the fragment domain | the back-shifted Fermi-gas systematics of von Egidy and Bucurescu, doi:10.1103/PhysRevC.72.044311, doi:10.1103/PhysRevC.73.049901, doi:10.1103/PhysRevC.80.054310, with the same mass evaluation |

Both are functions of those public evaluations alone, so they can be regenerated from them. Only
the header lines were rewritten when the naming convention was applied. `.gitattributes` keeps the
files exempt from line-ending normalisation.

## `wahl1988/`

Tables I–IV of Wahl (1988), doi:10.1016/0092-640X(88)90016-2; see the README there.

These are not inputs. They exist so that a change in the physics shows up as a failing test,
not as a number nobody checked.
