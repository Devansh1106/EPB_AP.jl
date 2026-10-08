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
