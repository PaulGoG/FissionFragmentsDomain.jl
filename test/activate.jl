using Pkg
Pkg.activate(@__DIR__; io = devnull)
# Resolve before instantiating: the package's own [compat] changes as it grows, and a manifest
# resolved against an older one would otherwise fail to load rather than update itself.
Pkg.resolve(; io = devnull)
Pkg.instantiate(; io = devnull)
