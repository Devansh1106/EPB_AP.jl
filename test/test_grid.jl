using EPB_AP: Grid, grad, divm, lap, logmean, harmonic, energy

@testset "grid and operators" begin
    g = Grid(0, 1, 37)
    @test g.dx ≈ 1 / 37
    @test g.x[1] ≈ g.dx / 2 && g.x[end] ≈ 1 - g.dx / 2
    p, q, v = rand(37), rand(37), rand(37)
    dx = g.dx
    # Proposition 2.2 (i) duality, (ii) self-adjointness, (iii) only constants in the kernel
    @test dx * sum(q .* divm(v, dx)) + dx * sum(grad(q, dx) .* v) ≈ 0 atol = 1e-12
    @test dx * sum(p .* lap(q, dx)) ≈ dx * sum(lap(p, dx) .* q)
    @test dx * sum(q .* lap(q, dx)) < 0
    @test maximum(abs, lap(fill(2.0, 37), dx)) == 0
    # (iv) second-order consistency of ΔM
    err(N) = (h = Grid(0, 1, N); maximum(abs, lap(sinpi.(2h.x), h.dx) .+ 4π^2 .* sinpi.(2h.x)))
    @test log2(err(50) / err(100)) ≈ 2 atol = 0.05
end

@testset "interface means" begin
    @test logmean(1.0, 2.0) ≈ 1 / log(2) rtol = 1e-15
    @test logmean(0.3, 0.3 + 1e-9) ≈ 0.3 + 5e-10 rtol = 1e-15
    @test logmean(1.0, 1.0 + 0.2) ≈ 0.2 / log(1.2) rtol = 1e-14   # just inside the series branch
    for (a, b) in [(1.0, 2.0), (0.3, 0.3 + 1e-9), (1e-4, 1.0), (2.0, 2.0)]
        @test logmean(a, b) * (log(b) - log(a)) ≈ b - a atol = 1e-15   # (10)
        @test min(a, b) <= logmean(a, b) <= max(a, b)
        @test 1 / harmonic(a, b) ≈ (1 / a + 1 / b) / 2
    end
end

@testset "discrete energy (14)" begin
    g = Grid(0, 1, 20)
    ρ, u, φ = fill(2.0, 20), fill(3.0, 20), fill(0.5, 20)
    @test energy(ρ, u, φ, 1.0, g.dx) ≈ 9 + exp(0.5) * (0.5 - 1)
end
