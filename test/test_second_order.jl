using EPB_AP: grad, avg, left, right, initial_potential, reconstruct, semidiscrete_rhs, faces, SecondOrder
using Random

@testset "minmod sign property (Proposition 4.1)" begin
    rng = MersenneTwister(4)
    w = randn(rng, 40)
    wm, wp, _ = reconstruct(w, 0.1, :minmod)
    β = (wp .- wm) ./ [w[right(i, 40)] - w[i] for i in 1:40]
    @test all(-1e-14 .<= β .<= 1 + 1e-14)
end

@testset "semidiscrete energy rate (93) equals (68) with the flux (61)" begin
    rng = MersenneTwister(20261008)
    for limiter in (:minmod, :none), c_D in (0.0, 1.0, 10.0), λ in (1.0, 1e-2)
        N, dx = 48, 1 / 48
        x = ((1:N) .- 0.5) .* dx
        ρ = 1 .+ 0.5sinpi.(2x) .+ (limiter == :minmod ? 0.3rand(rng, N) : 0.0)
        u = 0.3cospi.(2x) .+ (limiter == :minmod ? 0.5randn(rng, N) : 0.0)
        φ = initial_potential(ρ, λ, dx)
        sch = SecondOrder(; limiter, c_D)
        ρ̇, q̇, F, f = semidiscrete_rhs(ρ, ρ .* u, φ, dx, sch)
        rate = dx * sum(@. u * q̇ - u^2 / 2 * ρ̇ + φ * ρ̇)
        Δv = [u[right(i, N)] - u[i] for i in 1:N]
        σup = [F[e] >= 0 ? f.σv[e] : f.σv[right(e, N)] for e in 1:N]   # β of Remark 4.2
        pred = -sum(@. abs(F) * (Δv - dx * σup) * Δv) / 2 - dx * sum(grad(φ, dx) .* f.D)
        @test rate ≈ pred rtol = 1e-11
        @test pred <= 0
    end
end

@testset "mass flux is second-order consistent (Proposition 4.9)" begin
    for c_D in (0.0, 1.0)
        errs = map((50, 100, 200, 400)) do N
            dx = 1 / N
            x = ((1:N) .- 0.5) .* dx
            ρ, u = 1 .+ 0.3sinpi.(2x), 0.2cospi.(2x)
            ρ̇, = semidiscrete_rhs(ρ, ρ .* u, initial_potential(ρ, 1.0, dx), dx, SecondOrder(; limiter = :none, c_D))
            exact = @. -(0.6π * cospi(2x) * 0.2cospi(2x) - (1 + 0.3sinpi(2x)) * 0.4π * sinpi(2x))
            dx * sum(abs, ρ̇ .- exact)
        end
        @test all(1.95 .< log2.(errs[1:3] ./ errs[2:4]) .< 2.05)
    end
end

@testset "walls: no slope in the end cells, no flux through the walls" begin
    N, dx = 20, 0.05
    ρ, u = collect(range(1, 2; length = N)), collect(range(0.1, 0.5; length = N))
    _, _, σ = reconstruct(ρ, dx, :minmod; wall = true)
    @test σ[1] == σ[N] == 0
    ρ̇, q̇, F, f = semidiscrete_rhs(ρ, ρ .* u, zeros(N), dx, SecondOrder(); wall = true)
    @test F[end] == 0 && f.D[end] == 0
    @test abs(sum(ρ̇)) < 1e-12
end

@testset "SecondOrder options are checked" begin
    @test_throws ArgumentError SecondOrder(limiter = :vanleer)
end

@testset "fully discrete second order (§5.4 state, Δt = 0.4Δx)" begin
    run(N, λ, limiter) = solve(1 .+ 0.3sinpi.(2Grid(0, 1, N).x), 0.2cospi.(2Grid(0, 1, N).x), Grid(0, 1, N), λ, 0.2;
                               scheme = SecondOrder(; limiter), timestep = TimeStep(mode = :fixed, cfl = 0.4))
    restrict(v, N) = vec(sum(reshape(v, length(v) ÷ N, N); dims = 1)) ./ (length(v) ÷ N)
    # λ = 1 needs unlimited slopes: minmod clips the steep peak at x = 1/4 (Remark 4.10, order ≈ 1.5)
    for (λ, limiter) in ((1.0, :none), (1e-4, :minmod))
        ref = run(1600, λ, limiter)
        errs = map((50, 100, 200)) do N
            r = run(N, λ, limiter)
            @test r.status == :completed
            [sum(abs, f(r) .- restrict(f(ref), N)) / N for f in (r -> r.rho, r -> r.rho .* r.u, r -> r.phi)]
        end
        p = [log2.(errs[k] ./ errs[k + 1]) for k in 1:2]
        @test all(p[2] .> 1.8)
    end
end

@testset "second-order runs: mass, positivity, time-step rule" begin
    g = Grid(0, 1, 100)
    r = solve(1 .+ 0.3sinpi.(2g.x), 0.2cospi.(2g.x), g, 1.0, 0.3; scheme = SecondOrder())
    h = r.history
    @test r.status == :completed
    @test maximum(abs, h.mass .- h.mass[1]) <= 1e-13 * h.mass[1]
    @test all(h.min_rho .> 0)
    ρ, u = 1 .+ 0.3sinpi.(2g.x), 0.2cospi.(2g.x)
    @test 0 < EPB_AP.stable_dt(ρ, u, EPB_AP.initial_potential(ρ, 1.0, g.dx), g.dx, SecondOrder()) <= 0.4g.dx / 0.2   # (91)
    w = Grid(0, 1, 100; bc = :wall)
    r = solve(ifelse.(w.x .< 0.5, 1.0, 0.5), zero(w.x), w, 1e-2, 0.2; scheme = SecondOrder())
    @test r.status == :completed && abs(r.history.mass[end] - r.history.mass[1]) < 1e-13
end
