"""
    TimeStep(; mode = :cfl, cfl = 0.25, dtmax = Inf)

* `mode = :cfl`: Δt from the scheme's rule ([`stable_dt`](@ref)), checked a posteriori: first order with
  the computed flux, a violating step being repeated with Δt divided by the measured ratio; second order
  for positive densities, a failing step being repeated with Δt/2. Δt < 10⁻⁴Δx is a collapse.
* `mode = :fixed`: Δt = `cfl`·Δx, no check (§5.6, AP study); a second-order step with a non-positive
  density is a collapse.
* `dtmax`: upper bound on Δt/Δx (e.g. 0.4 in the order studies).
"""
Base.@kwdef struct TimeStep
    mode::Symbol = :cfl
    cfl::Float64 = 0.25
    dtmax::Float64 = Inf
end

"Result of [`solve`](@ref). `history` holds per-step time series; `snapshots` are `(t, ρ, u, φ)`."
struct Result
    grid::Grid
    rho::Vector{Float64}
    u::Vector{Float64}
    phi::Vector{Float64}
    t::Float64
    steps::Int
    status::Symbol
    history::NamedTuple
    snapshots::Vector{Any}
end

ap_residual(ρ, φ) = maximum(abs, exp.(φ) .- ρ)

# Scheme-specific parts of a step: the time-step rule, and one attempt returning `(ρ, u, φ, r)` where
# r > 1 asks for a repeat with Δt/r ((CFL)/(CFLθ) ratio for first order, 2 after a non-positive density
# for second order).
rule_dt(ρ, u, φ, dx, s::Scheme; wall) =
    stable_dt(ρ, u, φ, dx, eta_coefficient(ρ, s), kappa_coefficient(u, dx, s), s; wall)
rule_dt(ρ, u, φ, dx, s::SecondOrder; wall) = stable_dt(ρ, u, φ, dx, s; wall)

function attempt(ρ, u, φ, dt, dx, λ, s::Scheme; wall)
    ρ1, u1, φ1, F = step(ρ, u, φ, dt, dx, λ, eta_coefficient(ρ, s), kappa_coefficient(u, dx, s); wall)
    return ρ1, u1, φ1, cfl_ratio(ρ, F, dt, dx, s)
end

function attempt(ρ, u, φ, dt, dx, λ, s::SecondOrder; wall)
    try
        return step(ρ, u, φ, dt, dx, λ, s; wall)..., 1.0
    catch e
        e isa NonPositiveDensity || rethrow()
        return ρ, u, φ, 2.0
    end
end

"""
    solve(ρ0, u0, grid, λ, T; scheme = Scheme(), timestep = TimeStep(), snapshots = Float64[])

Integrate the scheme (`Scheme()` first order, `SecondOrder()` second order) from (ρ0, u0) to time T
on `grid` (periodic or walls);
φ⁰ is computed from ρ0 by (29).
`status` is `:completed` or `:collapse`.
"""
function solve(ρ0, u0, grid::Grid, λ, T; scheme = Scheme(), timestep = TimeStep(),
               snapshots = Float64[], maxretry = 100)
    dx, wall = grid.dx, grid.bc == :wall
    ρ, u = float.(copy(ρ0)), float.(copy(u0))
    φ = initial_potential(ρ, λ, dx; wall)
    hist = (t = [0.0], dt = [NaN], energy = [energy(ρ, u, φ, λ, dx; wall)], min_rho = [minimum(ρ)],
            mass = [dx * sum(ρ)], ap_residual = [ap_residual(ρ, φ)], cfl_ratio = [NaN])
    stops = sort(unique([filter(s -> 0 < s < T, snapshots); T]))
    snaps = Any[]
    t, n, status = 0.0, 0, :completed
    while t < T
        target = first(filter(s -> s > t, stops))
        dt = timestep.mode == :fixed ? timestep.cfl * dx : rule_dt(ρ, u, φ, dx, scheme; wall)
        dt = min(dt, timestep.dtmax * dx)
        if timestep.mode == :cfl && dt < 1e-4 * dx
            status = :collapse
            break
        end
        t + dt * (1 + 1e-6) >= target && (dt = target - t)   # land on target, no sliver step
        local ρ1, u1, φ1, r
        retried = false
        for _ in 1:maxretry
            ρ1, u1, φ1, r = attempt(ρ, u, φ, dt, dx, λ, scheme; wall)
            (timestep.mode == :fixed || r <= 1) && break
            dt /= r
            retried = true
        end
        failed = r > 1 && (timestep.mode == :cfl || scheme isa SecondOrder)   # fixed first order: no check
        if failed || (timestep.mode == :cfl && retried && dt < 1e-4 * dx)
            status = :collapse
            break
        end
        ρ, u, φ = ρ1, u1, φ1
        t = dt == target - t ? target : t + dt
        n += 1
        push!(hist.t, t); push!(hist.dt, dt); push!(hist.energy, energy(ρ, u, φ, λ, dx; wall))
        push!(hist.min_rho, minimum(ρ)); push!(hist.mass, dx * sum(ρ))
        push!(hist.ap_residual, ap_residual(ρ, φ)); push!(hist.cfl_ratio, r)
        t == target && target < T && push!(snaps, (t, copy(ρ), copy(u), copy(φ)))
    end
    return Result(grid, ρ, u, φ, t, n, status, hist, snaps)
end
