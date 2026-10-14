```@meta
CurrentModule = FissionFragmentsDomain
```

# Yields and averages

The fragment yield `Y(A, Z, TKE) = p(Z, A) Y(A, TKE)`, the probability measure on the
fragmentation domain, its one-dimensional marginal `Y(A)`, and the contractions that turn a
quantity tabulated over fragment configurations into an observable. A reconstruction or a
measured matrix may reach total kinetic energies the split cannot have; `bounded_fragment_yield`
cuts each `P(TKE|A,Z)` at `Q + E*_CN` and restores the measured moments over what remains.

```@autodocs
Modules = [FissionFragmentsDomain]
Pages = ["yields.jl", "bounded_yield.jl", "averaging.jl"]
```
