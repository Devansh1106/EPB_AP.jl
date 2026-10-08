using EPB_AP: Grid, grad, divm, lap, edgemap, avg, logmean, energy, right, left,
              solve_potential, initial_potential, step, eta_coefficient, kappa_coefficient,
              stable_dt, cfl_ratio, Scheme
using Random

# Random state of SPEC §9 item 8.
function random_state(rng, N)
    σ = 3rand(rng)
    lρ = σ .* randn(rng, N)
    rand(rng, Bool) && (lρ = (circshift(lρ, 1) .+ lρ .+ circshift(lρ, -1)) ./ 3)
    return exp.(lρ), 3rand(rng) .* randn(rng, N), 10.0^(5rand(rng) - 4)
end

@testset "Poisson solve (29) and (52)" begin
    g = Grid(0, 1, 64)
    ρ = 1 .+ 0.5sinpi.(2g.x)
    for λ in (1.0, 1e-3, 1e-6)
        φ = initial_potential(ρ, λ, g.dx)
        @test maximum(abs, -λ^2 .* lap(φ, g.dx) .+ exp.(φ) .- ρ) < 1e-10
    end
    μ, κ, λ, dt = rand(64), 0.3, 0.5, 0.01
    φ = solve_potential(zeros(64), ρ, μ, κ, λ, dt, g.dx)
    r = exp.(φ) .- λ^2 .* lap(φ, g.dx) .- dt .* divm(μ .* grad(φ, g.dx), g.dx) .+
        κ * λ^2 * dt .* lap(lap(φ, g.dx), g.dx) .- ρ
    @test maximum(abs, r) < 1e-8
end

@testset "semidiscrete energy rate: (93) equals (23) with D = 0" begin
    Random.seed!(20260923)
    g = Grid(0, 1, 200)
    x, dx, N = g.x, g.dx, 200
    states = [(1 .+ 0.3sinpi.(2x), 0.2cospi.(2x)), (1 .+ 0.8sinpi.(2x) .^ 3, 0.5cospi.(2x)),
              (1e-4 .+ (1 - 1e-4) .* exp.(-60 .* (x .- 0.5) .^ 2), zero(x)),
              (ifelse.(x .< 0.5, 1.0, 0.3), ifelse.(x .< 0.5, 0.4, -0.4)),
              (0.5 .+ rand(N), 2 .* rand(N) .- 1)]
    for (ρ, u) in states, λ in (1.0, 1e-2)
        φ = initial_potential(ρ, λ, dx)
        ρbar = edgemap(logmean, ρ)
        F = ρbar .* avg(u)
        G = [F[e] >= 0 ? u[e] * F[e] : u[right(e, N)] * F[e] for e in 1:N]
        a = ρbar .* grad(φ, dx)
        S = [-(a[i] + a[left(i, N)]) / 2 for i in 1:N]
        ρdot, mdot = -divm(F, dx), -divm(G, dx) .+ S
        rate93 = dx * sum(@. u * mdot - u^2 * ρdot / 2 + φ * ρdot)
        rate23 = -sum(abs.(F) .* (circshift(u, -1) .- u) .^ 2) / 2
        @test rate23 <= 0
        @test rate93 ≈ rate23 rtol = 1e-9 atol = 1e-12
    end
end

@testset "one step: conservation, (29), positivity, energy (Theorems 3.5, 3.16)" begin
    rng = Xoshiro(1)
    for kind in (:shift, :dissipation), _ in 1:200
        N = 48
        g = Grid(0, 1, N)
        ρ, u, λ = random_state(rng, N)
        θ = rand(rng, (0.2, 0.5, 0.8))
        sch = Scheme(; kind, theta = θ)
        φ = initial_potential(ρ, λ, g.dx)
        η = eta_coefficient(ρ, sch)
        κ = kind == :shift ? 0.0 : 10.0^(2rand(rng) - 1) * g.dx / 2 * (maximum(abs, u) + 1)
        dt = 1.0
        local ρ1, u1, φ1, F
        while true
            ρ1, u1, φ1, F = step(ρ, u, φ, dt, g.dx, λ, η, κ)
            cfl_ratio(ρ, F, dt, g.dx, sch) <= 1 && break
            dt /= 2
        end
        E0 = energy(ρ, u, φ, λ, g.dx)
        @test sum(ρ1) ≈ sum(ρ) rtol = 1e-12
        @test maximum(abs, -λ^2 .* lap(φ1, g.dx) .+ exp.(φ1) .- ρ1) <= 1e-8 * maximum(ρ1)
        @test all(ρ1 .>= (kind == :shift ? 4 / 5 : 1 - θ) .* ρ .* (1 - 1e-12))
        @test energy(ρ1, u1, φ1, λ, g.dx) <= E0 + 1e3 * eps() * max(1, abs(E0))
    end
end

@testset "time step rule satisfies the CFL condition" begin
    g = Grid(0, 1, 100)
    ρ, u = 1 .+ 0.5sinpi.(2g.x), cospi.(2g.x)
    for sch in (Scheme(), Scheme(kind = :dissipation))
        φ = initial_potential(ρ, 1.0, g.dx)
        η, κ = eta_coefficient(ρ, sch), kappa_coefficient(u, g.dx, sch)
        dt = stable_dt(ρ, u, φ, g.dx, η, κ, sch)
        @test 0 < dt < Inf
        F = step(ρ, u, φ, dt, g.dx, 1.0, η, κ)[4]
        @test cfl_ratio(ρ, F, dt, g.dx, sch) <= 1.05
    end
end

@testset "scheme kind is checked" begin
    @test Scheme(kind = :dissipation).kind == :dissipation
    @test_throws ArgumentError Scheme(kind = :dissipative)
end
