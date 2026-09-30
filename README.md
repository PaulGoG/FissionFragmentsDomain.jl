# FissionFragmentsDomain.jl

[![CI](https://github.com/PaulGoG/FissionFragmentsDomain.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/PaulGoG/FissionFragmentsDomain.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation (stable)](https://img.shields.io/badge/docs-stable-blue.svg)](https://PaulGoG.github.io/FissionFragmentsDomain.jl/stable/)
[![Documentation (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://PaulGoG.github.io/FissionFragmentsDomain.jl/dev/)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

```
FissionFragmentsDomain.jl/
├── src/            the library: systems, masses, charge, domain, level density, energetics,
│                   temperature ratio, yields, averages, curves and run records
├── data/           the shipped AME2020 mass table and the formats every reader accepts
├── test/           unit tests, reference values from the archived ²⁵²Cf run, static QA
├── docs/           Documenter site
├── activate.jl     activates and instantiates the package environment
└── check.jl        formatter, then the test suite
```

The fragments of binary fission as models of prompt emission and of the excitation-energy
partition see them. A fissioning system splits into complementary pairs
`(A_L, Z_L) + (A_H, Z_H)` at a total kinetic energy `TKE`; this package holds what such models
share, so that a calculation extracting the temperature ratio `R_T(A_H)` from measured `ν(A)` and
one applying it to the emission of prompt neutrons work on one domain rather than two copies that
drift apart:

- nuclides, fissioning systems and the system token `Cf252_sf` that names them everywhere;
- the 2020 atomic mass evaluation, shipped, with `Q`-values and separation energies;
- the isobaric charge distribution — evaluated tables, or the Wahl `Zₚ` model with per-reaction
  or systematic parameters — and the fragmentation domain built from it;
- fragment level density parameters (back-shifted Fermi gas, Gilbert–Cameron), their ratio across
  a pair averaged over the charge distribution, and the relation `E*_H/TXE = 1/(1 + R_a R_T²)` with
  its inverse;
- the energy balance `TXE = Q + E*_CN − TKE`;
- the fragment yield `Y(A, Z, TKE) = p(Z, A) Y(A, TKE)`, its marginal `Y(A)`, and averages over
  it;
- curves of the heavy-fragment mass, `R_T(A_H)` and `E*_H/TXE (A_H)`, and the run record through
  which an extraction hands them on.

It is consumed by [`DeterministicSequentialEmission.jl`](https://github.com/PaulGoG/DeterministicSequentialEmission.jl)
and [`FissionTemperatureRatio.jl`](https://github.com/PaulGoG/FissionTemperatureRatio.jl).

## Getting started

Julia 1.11 or later. The package is installed from its repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/PaulGoG/FissionFragmentsDomain.jl")
```

A consuming package declares it under `[sources]` in its `Project.toml`:

```toml
[sources]
FissionFragmentsDomain = {url = "https://github.com/PaulGoG/FissionFragmentsDomain.jl"}
```

From a clone, every environment activates and instantiates itself:

```bash
julia -i activate.jl               # REPL in the package environment
julia check.jl                     # format, then run the tests
julia docs/make.jl                 # build the documentation into docs/build
```

Every test rests on a public source shipped with the package, except two that need the
Gilbert–Cameron shell-correction table. Those read a data directory laid out as `data/README.md`
describes, and are reported as skipped where it is absent:

```bash
FISSION_FRAGMENTS_DOMAIN_TEST_DATA=/path/to/data julia check.jl
```

```julia
using FissionFragmentsDomain

system = spontaneous_fission(Nuclide(98, 252))
masses = read_mass_excess_table(AME2020_MASS_EXCESS_FILE)
charge = charge_model(masses, system)                           # Wahl (1988), CF252S
domain = fragmentation_domain(system, charge, 126:174)
R_a = level_density_ratio(RatioOfMeans(), BackShiftedFermiGas(masses), domain)
heavy_excitation_fraction(1.2, R_a[140])                        # E*_H/TXE at A_H = 140
```

## Status

| Component | State |
|---|---|
| Systems, masses, `Q` | complete; AME2020 shipped |
| Charge distribution | evaluated tables; Wahl (1988) for its four reactions, reproducing its Tables I–IV row by row; Wahl (2002) systematics, low-energy branch; conventional means |
| Fragmentation domain | complete; optional zero polarization at the symmetric split |
| Level density | back-shifted Fermi gas (von Egidy–Bucurescu), Gilbert–Cameron with its deformed branch |
| Temperature ratio | forward and inverse relation, `R_a` in either averaging order |
| Yields | joint `Y(A, TKE)`, factorized reconstruction, `Y(A)` |
| Run records | reader; the writer follows with the rebuilt producer |

<details>
<summary>Full tree</summary>

```
FissionFragmentsDomain.jl/
├── .github/                 CI workflow and Dependabot configuration
├── data/
│   ├── README.md             reader formats, provenance, lineage of the charge tables
│   └── reference/mass_excess_ame2020.dat
├── docs/
│   ├── make.jl, activate.jl, Project.toml
│   └── src/                  index, system, fragmentation, energetics, yields, curves,
│                             references, refs.bib
├── formatter/                pinned JuliaFormatter environment
├── src/
│   ├── FissionFragmentsDomain.jl
│   ├── nuclides.jl           nuclides, systems, the system token and record
│   ├── tables.jl             positional table reader
│   ├── masses.jl             mass excesses, separation energies, Q
│   ├── charge.jl             isobaric charge distribution
│   ├── wahl.jl               Wahl (2002) Zₚ systematics
│   ├── zp.jl                 the Zₚ model: eq. (7), products and fragments
│   ├── wahl1988.jl           Wahl (1988) per-reaction model, layer resolution
│   ├── domain.jl             fragmentation domain, splits, charge windows
│   ├── leveldensity.jl       BSFG and Gilbert–Cameron
│   ├── energetics.jl         E*_CN, TXE, kinematics, energetics sweep
│   ├── temperature_ratio.jl  R_a, R_T ↔ E*_H/TXE
│   ├── yields.jl             Y(A,TKE), factorized yield, Y(A)
│   ├── averaging.jl          yield-weighted contractions
│   ├── curves.jl             curves of A_H
│   └── manifest.jl           temperature-ratio run record
├── test/
│   ├── runtests.jl, physics.jl, activate.jl, Project.toml
│   ├── test_*.jl
│   └── references/           Q and a over the ²⁵²Cf domain; Wahl (1988), Tables I–IV
├── activate.jl
├── check.jl
├── CHANGELOG.md
├── CITATION.cff
└── LICENSE
```

</details>

## How to cite

The software is citable through `CITATION.cff`:

```bibtex
@software{Gogita_FissionFragmentsDomain_jl,
  author  = {Gogîță, Paul-Adrian},
  title   = {FissionFragmentsDomain.jl},
  url     = {https://github.com/PaulGoG/FissionFragmentsDomain.jl},
  license = {MIT},
}
```

The mass evaluation, the charge distribution systematics and the level density systematics each
carry their own citation; they are listed with DOIs in the documentation's References page and in
`data/README.md`.
