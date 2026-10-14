"""
    kinetic_energy_limit(table, system, fragment) -> Union{Float64,Nothing}

The physical upper limit of the total kinetic energy of the split that produces `fragment`, in
MeV: `Q + E*_CN`, where `TXE = Q + E*_CN − TKE` vanishes. A configuration at or beyond it has no
excitation energy to share and is not a fission event. `nothing` where the mass table cannot form
`Q` or `E*_CN`.
"""
function kinetic_energy_limit(
    table::MassExcessTable,
    system::FissioningSystem,
    fragment::Nuclide,
)
    Q = q_value(table, system, fragment)
    Q === nothing && return nothing
    E_compound = compound_nucleus_excitation(table, system)
    E_compound === nothing && return nothing
    return value(Q + E_compound)
end

"""
    maximum_kinetic_energy_limit(table, system, domain) -> NamedTuple

The largest [`kinetic_energy_limit`](@ref) over the fragmentations of `domain`, as
`(limit, fragmentation)`: the total kinetic energy a grid must reach so that it cuts no physical
part of any `P(TKE|A,Z)` the domain spans. Throws an `ArgumentError` where no fragmentation of the
domain has a limit.
"""
function maximum_kinetic_energy_limit(
    table::MassExcessTable,
    system::FissioningSystem,
    domain::FragmentationDomain,
)
    best = -Inf
    where_ = nothing
    for fragmentation in domain
        limit = kinetic_energy_limit(table, system, fragmentation.heavy)
        limit === nothing && continue
        if limit > best
            best = limit
            where_ = fragmentation
        end
    end
    where_ === nothing && throw(
        ArgumentError("no fragmentation of the domain has a Q-value in $(table.source)"),
    )
    return (limit = best, fragmentation = where_)
end

"""
    KINETIC_ENERGY_MOMENTS

The moments of `P(TKE|A)` that [`bounded_fragment_yield`](@ref) restores after bounding it at the
physical limit: `:none` keeps the bounded distribution as it falls, `:mean` restores the measured
`⟨TKE⟩(A)`, `:mean_and_width` restores `⟨TKE⟩(A)` and `σ_TKE(A)` together.
"""
const KINETIC_ENERGY_MOMENTS = (:none, :mean, :mean_and_width)

"""
    KineticEnergyBound

How the kinetic energy distribution of one fragment mass came out of
[`bounded_fragment_yield`](@ref), all energies in MeV:

- `A` — the fragment mass;
- `measured_mean`, `measured_dispersion` — the `⟨TKE⟩(A)` and `σ_TKE(A)` the bounded distribution
  was asked to keep (`NaN` where no width was measured);
- `mean`, `dispersion` — the mean and standard deviation over the retained cells, charges pooled;
- `parent_mean`, `parent_dispersion` — the Gaussian the retained cells are cut from, `μ′` and `σ′`,
  where the row is a factorized Gaussian (`NaN` for a measured joint row or a mass placed at its
  mean; `Inf` where no Gaussian has the measured width once bounded);
- `removed` — the share of the unbounded row, charge-weighted, at or beyond `Q + E*_CN`.
"""
struct KineticEnergyBound
    A::Int
    measured_mean::Float64
    measured_dispersion::Float64
    mean::Float64
    dispersion::Float64
    parent_mean::Float64
    parent_dispersion::Float64
    removed::Float64
end

function Base.show(io::IO, bound::KineticEnergyBound)
    print(
        io,
        "KineticEnergyBound(A=",
        bound.A,
        ", ⟨TKE⟩ ",
        round(bound.mean; digits = 4),
        " against ",
        round(bound.measured_mean; digits = 4),
        " MeV, ",
        round(100 * bound.removed; sigdigits = 3),
        " % removed)",
    )
    return nothing
end

"""
    bounded_fragment_yield(domain, distribution, table, experimental; moments = :mean,
                           measured = nothing, total = DEFAULT_YIELD_TOTAL)
        -> (FragmentYield, Vector{KineticEnergyBound})

Spread a `Y(A, TKE)` over charge numbers as [`build_fragment_yield`](@ref) does, but bounded at
the physical limit of each charge split:

```
Y(A, Z, TKE) = Y(A) p(Z, A) P(TKE | A, Z),    P(TKE | A, Z) = 0 for TKE ≥ Q(A, Z) + E*_CN
```

A configuration at or beyond [`kinetic_energy_limit`](@ref) is impossible, not merely
unmeasured, so its weight belongs to the retained configurations of the same `(A, Z)`. Each
`P(TKE|A,Z)` is the mass row of `experimental` cut at the limit and normalised over what remains,
so the mass yield `Y(A)` and the charge distribution `p(Z|A)`, both integrals over `TKE`, keep
their measured values. A bound alone lowers the mean kinetic energy of a mass whose row reaches
the limit; `moments` restores what was measured:

- `:none` — the bounded rows as they fall;
- `:mean` — every row of mass `A` multiplied by `exp(λ (TKE − ⟨TKE⟩(A)))`, one `λ` per mass,
  solved so that the mean over the retained cells, charges pooled, equals `⟨TKE⟩(A)`;
- `:mean_and_width` — the factor `exp(λ₁ x + λ₂ x²)`, `x = TKE − ⟨TKE⟩(A)`, with `λ₁` and `λ₂`
  solved so that the standard deviation equals `σ_TKE(A)` as well, where one was measured.

Of all distributions with the measured `Y(A)`, `p(Z|A)` and moments, this is the one closest to
the bounded row in the Kullback–Leibler sense: the minimum-discrimination-information solution
under linear constraints is an exponential tilt (Kullback, *Information Theory and Statistics*,
Wiley (1959); Csiszár, *Ann. Probab.* **3**, 146 (1975), doi:10.1214/aop/1176996454). For a
Gaussian row the tilt moves the mean and the width of the Gaussian itself, which is how the
method on the marginal files reports it. The mean of a mass rises monotonically with `λ`, so each
solve is a one-dimensional bracketed root; the width, at a fixed mean, rises monotonically with
`λ₂`.

`measured` holds the `(⟨TKE⟩(A), σ_TKE(A))` to restore, keyed by mass, a width of `nothing` where
none was measured; `nothing` takes each row's own moments, which is what a measured joint
distribution reports. The result is renormalized to `total`, dividing out once the charge weight
the retained charges and the mass table do not reach, as [`build_fragment_yield`](@ref) does. A
mass whose measured mean the retained cells cannot reach throws an `ArgumentError`.
"""
function bounded_fragment_yield(
    domain::FragmentationDomain,
    distribution::ChargeModel,
    table::MassExcessTable,
    experimental::MassEnergyYield;
    moments::Symbol = :mean,
    measured::Union{Nothing, AbstractDict} = nothing,
    total::Real = DEFAULT_YIELD_TOTAL,
)
    values, bounds, _ = _bounded_fragment_yield(
        domain,
        distribution,
        table,
        experimental,
        moments,
        measured,
        total,
    )
    return FragmentYield(values, Float64(total)), bounds
end

"""
    bounded_fragment_yield(domain, distribution, table, mass_yield_path, kinetic_energy_path,
                           energies; dispersion = nothing, symmetrize = false, moments = :mean,
                           total = DEFAULT_YIELD_TOTAL)
        -> (FragmentYield, Vector{KineticEnergyBound})

The bounded yield of a [`factorized_yield`](@ref) reconstruction, from the same marginal files
and keywords: each mass row the Gaussian of the measured `⟨TKE⟩(A)` and `σ_TKE(A)` on `energies`,
bounded at `Q + E*_CN` per charge, and re-centred (`moments = :mean`, the width kept) or
re-centred and re-scaled (`:mean_and_width`) so that the bounded distribution of each mass keeps
the measured moments; see the method on a [`MassEnergyYield`](@ref). The records carry the
Gaussian the retained cells are cut from, `μ′` and `σ′`.

The grid is part of the reconstruction here, not of a measurement, so it must not truncate a
distribution where it is physical: an `ArgumentError` is thrown where `energies`, continued by
its last step, stops short of the limit of a fragment that carries yield. The lower end is not
checked; the restored moments are those over the grid.
"""
function bounded_fragment_yield(
    domain::FragmentationDomain,
    distribution::ChargeModel,
    table::MassExcessTable,
    mass_yield_path::AbstractString,
    kinetic_energy_path::AbstractString,
    energies;
    dispersion::Union{Nothing, Real, KineticEnergyDispersion} = nothing,
    symmetrize::Bool = false,
    moments::Symbol = :mean,
    total::Real = DEFAULT_YIELD_TOTAL,
)
    compound_mass = domain.system.compound.A
    rows = factorized_yield(
        mass_yield_path,
        kinetic_energy_path,
        energies,
        compound_mass;
        dispersion = dispersion,
        total = total,
        symmetrize = symmetrize,
    )
    measured = kinetic_energy_moments(
        mass_yield_path,
        kinetic_energy_path,
        compound_mass;
        dispersion = dispersion,
        symmetrize = symmetrize,
    )
    _check_grid_reaches_limit(domain, distribution, table, rows)
    values, bounds, tilts =
        _bounded_fragment_yield(domain, distribution, table, rows, moments, measured, total)
    # A factorized row is the Gaussian of the measured moments sampled on the grid, so the tilt
    # exp(λ₁x + λ₂x²), x = TKE − μ, is that Gaussian moved and rescaled:
    # 1/σ′² = 1/σ² − 2λ₂, μ′ = μ + λ₁σ′².
    for (index, bound) in pairs(bounds)
        μ, σ = bound.measured_mean, bound.measured_dispersion
        isfinite(σ) || continue
        λ₁, λ₂ = tilts[index]
        curvature = 1 / σ^2 - 2λ₂
        σ′ = curvature > 0 ? sqrt(1 / curvature) : Inf
        μ′ = isfinite(σ′) ? μ + λ₁ * σ′^2 : Inf
        bounds[index] = KineticEnergyBound(
            bound.A,
            bound.measured_mean,
            bound.measured_dispersion,
            bound.mean,
            bound.dispersion,
            μ′,
            σ′,
            bound.removed,
        )
    end
    return FragmentYield(values, Float64(total)), bounds
end

# A reconstruction grid that ends below the limit of a fragment carrying yield cuts that
# fragment's distribution where it is physical.
function _check_grid_reaches_limit(
    domain::FragmentationDomain,
    distribution::ChargeModel,
    table::MassExcessTable,
    rows::MassEnergyYield,
)
    grid = rows.energies
    last_step = length(grid) >= 2 ? grid[end] - grid[end - 1] : 0.0
    reach = grid[end] + last_step
    carrying = Set(rows.masses)
    worst = -Inf
    culprit = nothing
    for nuclide in keys(fragment_charge_probabilities(domain, distribution))
        nuclide.A in carrying || continue
        limit = kinetic_energy_limit(table, domain.system, nuclide)
        limit === nothing && continue
        if limit > reach && limit > worst
            worst = limit
            culprit = nuclide
        end
    end
    culprit === nothing && return nothing
    throw(
        ArgumentError(
            "the kinetic energy grid ends at $(grid[end]) MeV, and its next step would reach \
             $reach MeV, short of the physical limit Q + E*_CN = $(round(worst; digits = 3)) \
             MeV of $culprit; a bounded reconstruction needs the grid to reach the limit, to \
             at least $(round(worst - last_step; digits = 3)) MeV",
        ),
    )
end

# One charge of one mass row: its probability and the grid indices of the cells it keeps.
const _RetainedCharge = Tuple{Nuclide, Float64, Vector{Int}}

function _bounded_fragment_yield(
    domain::FragmentationDomain,
    distribution::ChargeModel,
    table::MassExcessTable,
    experimental::MassEnergyYield,
    moments::Symbol,
    measured::Union{Nothing, AbstractDict},
    total::Real,
)
    moments in KINETIC_ENERGY_MOMENTS || throw(
        ArgumentError(
            "moments must be one of $(join(map(repr, KINETIC_ENERGY_MOMENTS), ", ")), got \
             $(repr(moments))",
        ),
    )
    total > 0 ||
        throw(ArgumentError("the yield normalization must be positive, got $total"))
    system = domain.system
    E_compound = compound_nucleus_excitation(table, system)
    E_compound === nothing && throw(
        ArgumentError(
            "the neutron separation energy of the compound nucleus $(system.compound) is not \
             available from $(table.source); the excitation energy cannot be formed",
        ),
    )

    by_mass = Dict{Int, Vector{Tuple{Nuclide, Float64}}}()
    for (nuclide, p) in fragment_charge_probabilities(domain, distribution)
        push!(get!(by_mass, nuclide.A, Tuple{Nuclide, Float64}[]), (nuclide, p))
    end

    grid = experimental.energies
    values = Dict{Tuple{Int, Int, Float64}, Float64}()
    bounds = KineticEnergyBound[]
    tilts = Tuple{Float64, Float64}[]
    running = 0.0
    for A in sort(collect(keys(by_mass)))
        present = [haskey(experimental.values, (A, TKE)) for TKE in grid]
        any(present) || continue
        q = [present[j] ? experimental.values[(A, grid[j])] : 0.0 for j in eachindex(grid)]
        row_total = sum(q)
        row_total > 0 || continue

        target_mean, target_dispersion = _row_targets(grid, q, present, measured, A)
        target_mean === nothing && continue

        charges = _RetainedCharge[]
        full = 0.0
        kept = 0.0
        for (nuclide, p) in sort(by_mass[A]; by = entry -> entry[1].Z)
            Q = q_value(table, system, nuclide)
            Q === nothing && continue
            full += p * row_total
            retained = Int[
                j for j in eachindex(grid) if present[j] &&
                    total_excitation_energy(Q, E_compound, grid[j]) !== nothing
            ]
            isempty(retained) && continue
            weight = sum(q[j] for j in retained)
            weight > 0 || continue
            kept += p * weight
            push!(charges, (nuclide, p, retained))
        end
        isempty(charges) && continue

        x = [TKE - target_mean for TKE in grid]
        ℓ = [q[j] > 0 ? log(q[j]) : -Inf for j in eachindex(grid)]
        λ₁, λ₂ = _restore_moments(x, ℓ, charges, moments, target_dispersion, A)
        M, V, _ = _tilted_moments(x, ℓ, charges, λ₁, λ₂)

        for (nuclide, p, retained) in charges
            exponents = [ℓ[j] + λ₁ * x[j] + λ₂ * x[j]^2 for j in retained]
            shift = maximum(exponents)
            weights = exp.(exponents .- shift)
            norm = sum(weights)
            for (k, j) in pairs(retained)
                cell = row_total * p * weights[k] / norm
                values[(A, nuclide.Z, grid[j])] = cell
                running += cell
            end
        end
        push!(
            bounds,
            KineticEnergyBound(
                A,
                target_mean,
                target_dispersion === nothing ? NaN : target_dispersion,
                target_mean + M,
                sqrt(max(V, 0.0)),
                NaN,
                NaN,
                full > 0 ? 1 - kept / full : 0.0,
            ),
        )
        push!(tilts, (λ₁, λ₂))
    end
    isempty(values) && throw(
        ArgumentError("no fragment of the domain keeps a yield below its physical limit"),
    )
    scale = total / running
    for key in keys(values)
        values[key] *= scale
    end
    return values, bounds, tilts
end

# The moments a row is to keep: the measured ones where given, else the row's own over every
# cell it holds, before the bound.
function _row_targets(grid, q, present, measured, A::Integer)
    if measured === nothing
        weight = sum(q)
        mean = sum(grid[j] * q[j] for j in eachindex(grid) if present[j]) / weight
        variance =
            sum((grid[j] - mean)^2 * q[j] for j in eachindex(grid) if present[j]) / weight
        return mean, sqrt(variance)
    end
    entry = get(measured, A, nothing)
    entry === nothing && return nothing, nothing
    mean, dispersion = entry
    return Float64(mean), dispersion === nothing ? nothing : Float64(dispersion)
end

# Mean and variance of x over the retained cells of a mass, charges pooled with their weights
# p, each charge normalised on its own cells, under the tilt exp(λ₁x + λ₂x²); and dM/dλ₁, the
# charge-weighted mean of the within-charge variances.
function _tilted_moments(
    x::Vector{Float64},
    ℓ::Vector{Float64},
    charges::Vector{_RetainedCharge},
    λ₁::Float64,
    λ₂::Float64,
)
    P = 0.0
    first_moment = 0.0
    second_moment = 0.0
    slope = 0.0
    for (_, p, retained) in charges
        shift = -Inf
        for j in retained
            exponent = ℓ[j] + λ₁ * x[j] + λ₂ * x[j]^2
            exponent > shift && (shift = exponent)
        end
        isfinite(shift) || continue
        norm = 0.0
        m₁ = 0.0
        m₂ = 0.0
        for j in retained
            w = exp(ℓ[j] + λ₁ * x[j] + λ₂ * x[j]^2 - shift)
            norm += w
            m₁ += w * x[j]
            m₂ += w * x[j]^2
        end
        m₁ /= norm
        m₂ /= norm
        P += p
        first_moment += p * m₁
        second_moment += p * m₂
        slope += p * (m₂ - m₁^2)
    end
    M = first_moment / P
    return M, second_moment / P - M^2, slope / P
end

function _restore_moments(x, ℓ, charges, moments::Symbol, target_dispersion, A::Integer)
    moments === :none && return 0.0, 0.0
    if moments === :mean || target_dispersion === nothing
        return _solve_location(x, ℓ, charges, 0.0, A), 0.0
    end
    variance = Float64(target_dispersion)^2
    excess(λ₂) =
        _tilted_moments(x, ℓ, charges, _solve_location(x, ℓ, charges, λ₂, A), λ₂)[2] -
        variance
    # The width at the measured mean rises with λ₂; bracket the root by doubling, then bisect.
    at_zero = excess(0.0)
    at_zero == 0 && return _solve_location(x, ℓ, charges, 0.0, A), 0.0
    step = 0.05 / variance
    if at_zero < 0
        lower, upper = 0.0, step
        for _ in 1:80
            excess(upper) > 0 && break
            lower, upper = upper, 2upper
        end
    else
        lower, upper = -step, 0.0
        for _ in 1:80
            excess(lower) < 0 && break
            lower, upper = 2lower, lower
        end
    end
    excess(lower) <= 0 <= excess(upper) || throw(
        ArgumentError(
            "the bounded kinetic energy distribution of A = $A cannot take the measured \
             width $(target_dispersion) MeV at its measured mean",
        ),
    )
    for _ in 1:200
        middle = (lower + upper) / 2
        (middle == lower || middle == upper) && break
        gap = excess(middle)
        if abs(gap) <= 1e-12 * variance
            lower = upper = middle
            break
        end
        gap < 0 ? (lower = middle) : (upper = middle)
    end
    λ₂ = (lower + upper) / 2
    return _solve_location(x, ℓ, charges, λ₂, A), λ₂
end

# λ₁ at fixed λ₂ such that the pooled mean of x vanishes: the mean rises monotonically with λ₁
# (its derivative is the pooled within-charge variance), so a bracket is found by doubling and
# refined by Newton steps that fall back to bisection.
function _solve_location(x, ℓ, charges, λ₂::Float64, A::Integer)
    M, _, slope = _tilted_moments(x, ℓ, charges, 0.0, λ₂)
    abs(M) <= 1e-12 && return 0.0
    slope > 0 || throw(
        ArgumentError(
            "the retained cells of A = $A cannot move their mean; the measured ⟨TKE⟩(A) is \
             out of reach",
        ),
    )
    step = abs(M) / slope
    lower, upper = M > 0 ? (-step, 0.0) : (0.0, step)
    for _ in 1:200
        M > 0 ? (_tilted_moments(x, ℓ, charges, lower, λ₂)[1] < 0 && break) :
        (_tilted_moments(x, ℓ, charges, upper, λ₂)[1] > 0 && break)
        M > 0 ? (lower *= 2) : (upper *= 2)
    end
    (
        _tilted_moments(x, ℓ, charges, lower, λ₂)[1] <=
        0 <=
        _tilted_moments(x, ℓ, charges, upper, λ₂)[1]
    ) || throw(
        ArgumentError(
            "the measured ⟨TKE⟩(A) of A = $A lies outside the reach of the cells below the \
             physical limit",
        ),
    )
    λ = -M / slope
    (lower < λ < upper) || (λ = (lower + upper) / 2)
    for _ in 1:200
        M, _, slope = _tilted_moments(x, ℓ, charges, λ, λ₂)
        abs(M) <= 1e-12 && return λ
        M < 0 ? (lower = λ) : (upper = λ)
        candidate = slope > 0 ? λ - M / slope : (lower + upper) / 2
        λ = lower < candidate < upper ? candidate : (lower + upper) / 2
        (λ == lower || λ == upper) && return λ
    end
    return λ
end
