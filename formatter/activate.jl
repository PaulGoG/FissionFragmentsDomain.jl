# Activate and instantiate the formatter environment, silently.
#
# Kept apart from the package and from the test environment: the formatter is a tool the code is
# checked with, not a dependency of anything, and pinning it here means the gate CI applies and
# the gate `check.jl` applies are the same version.
using Pkg
Pkg.activate(@__DIR__; io = devnull)
Pkg.resolve(; io = devnull)
Pkg.instantiate(; io = devnull)
