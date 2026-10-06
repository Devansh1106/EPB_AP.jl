"""
    TimeStep(; mode = :cfl, cfl = 0.25, dtmax = Inf)

* `mode = :cfl`: Δt from the scheme's rule ([`stable_dt`](@ref)), checked a posteriori with the computed
  flux; a violating step is repeated with Δt divided by the measured ratio. Δt < 10⁻⁴Δx is a collapse.
* `mode = :fixed`: Δt = `cfl`·Δx, no check (§5.6, AP study).
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

"""
    solve(ρ0, u0, grid, ε, T; scheme = Scheme(), timestep = TimeStep(), snapshots = Float64[])

Integrate the first-order scheme from (ρ0, u0) to time T; φ⁰ is computed from ρ0 by (29).
`status` is `:completed` or `:collapse`.
"""
function solve(ρ0, u0, grid::Grid, ε, T; scheme = Scheme(), timestep = TimeStep(),
               snapshots = Float64[], maxretry = 100)
    dx = grid.dx
    ρ, u = float.(copy(ρ0)), float.(copy(u0))
    φ = initial_potential(ρ, ε, dx)
    hist = (t = [0.0], dt = [NaN], energy = [energy(ρ, u, φ, ε, dx)], min_rho = [minimum(ρ)],
            mass = [dx * sum(ρ)], ap_residual = [ap_residual(ρ, φ)], cfl_ratio = [NaN])
    stops = sort(unique([filter(s -> 0 < s < T, snapshots); T]))
    snaps = Any[]
    t, n, status = 0.0, 0, :completed
    while t < T
        target = first(filter(s -> s > t, stops))
        η = eta_coefficient(ρ, scheme)
        κ = kappa_coefficient(u, dx, scheme)
        dt = timestep.mode == :fixed ? timestep.cfl * dx : stable_dt(ρ, u, φ, dx, η, κ, scheme)
        dt = min(dt, timestep.dtmax * dx)
        if timestep.mode == :cfl && dt < 1e-4 * dx
            status = :collapse
            break
        end
        t + dt * (1 + 1e-6) >= target && (dt = target - t)   # land on target, no sliver step
        local ρ1, u1, φ1, r
        retried = false
        for _ in 1:maxretry
            ρ1, u1, φ1, F = step(ρ, u, φ, dt, dx, ε, η, κ)
            r = cfl_ratio(ρ, F, dt, dx, scheme)
            (timestep.mode == :fixed || r <= 1) && break
            dt /= r
            retried = true
        end
        if timestep.mode == :cfl && (r > 1 || (retried && dt < 1e-4 * dx))
            status = :collapse
            break
        end
        ρ, u, φ = ρ1, u1, φ1
        t = dt == target - t ? target : t + dt
        n += 1
        push!(hist.t, t); push!(hist.dt, dt); push!(hist.energy, energy(ρ, u, φ, ε, dx))
        push!(hist.min_rho, minimum(ρ)); push!(hist.mass, dx * sum(ρ))
        push!(hist.ap_residual, ap_residual(ρ, φ)); push!(hist.cfl_ratio, r)
        t == target && target < T && push!(snaps, (t, copy(ρ), copy(u), copy(φ)))
    end
    return Result(grid, ρ, u, φ, t, n, status, hist, snaps)
end
