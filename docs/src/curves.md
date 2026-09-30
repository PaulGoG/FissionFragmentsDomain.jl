```@meta
CurrentModule = FissionFragmentsDomain
```

# Curves and run records

Functions of the heavy-fragment mass, `R_T(A_H)` and `E*_H/TXE (A_H)`, and the run record through
which a temperature-ratio extraction hands its curves to a consumer. The record declares what each
file tabulates, so that `R_T` and the multiplicity ratio `r_ν`, which occupy the same column of two
files written side by side, cannot be confused.

```@autodocs
Modules = [FissionFragmentsDomain]
Pages = ["curves.jl", "manifest.jl"]
```
