using EPB_AP: CASES, lap, riemann_speeds, SolitonProfile

@testset "test case definitions" begin
    for name in CASES
        c = testcase(name)
        @test c.name == name
        g = Grid(c.domain..., 64)
        ρ, u = c.init(g.x)
        @test length(ρ) == length(u) == 64 && all(ρ .> 0)
    end
    @test testcase(:smooth; eps = 1e-2, N = 50).eps == 1e-2
    @test_throws ArgumentError testcase(:nope)
end

@testset "exact solutions" begin
    um, us = riemann_speeds(0.5)
    @test um ≈ 0.345708 atol = 1e-6               # [1, (7.9)], §5.5
    @test us ≈ 1.17786117 atol = 1e-7
    # soliton: Poisson residual of the profile is the O(Δx²) error of ΔM, and it is a travelling wave
    c = testcase(:soliton)
    res(N) = (g = Grid(0, 40, N); (ρ, u, φ) = c.exact(g.x, 0.0);
              maximum(abs, -lap(φ, g.dx) .+ exp.(φ) .- ρ))
    @test log2(res(1000) / res(2000)) ≈ 2 atol = 0.1
    g = Grid(0, 40, 400)
    @test c.exact(g.x, 40 / 1.3)[3] ≈ c.exact(g.x, 0.0)[3] atol = 1e-12
    # simple wave: constant along ẋ = u + 1
    ρ, u, _ = EPB_AP.simplewave_exact([0.3], 0.0)
    @test EPB_AP.simplewave_exact([0.3 + 0.2 * (u[1] + 1)], 0.2)[2] ≈ u atol = 1e-13
end

@testset "properties along runs: mass, positivity, energy" begin
    for (name, kind) in [(:smooth, :shift), (:steep, :shift), (:nearvacuum, :shift),
                         (:bump, :dissipation), (:step, :dissipation), (:expansion, :dissipation)]
        r = solve(testcase(name; T = 0.2, scheme = Scheme(kind = kind)))
        h = r.history
        @test r.status == :completed
        @test maximum(abs, h.mass .- h.mass[1]) <= 1e-13 * h.mass[1]
        @test all(h.min_rho .> 0)
        @test all(diff(h.energy) .<= 1e3 * eps() * max(1, abs(h.energy[1])))
        @test all(filter(!isnan, h.cfl_ratio) .<= 1)
    end
end

@testset "asymptotic preservation: residual ∝ ε²" begin
    for (name, kind) in [(:smooth, :shift), (:ap, :dissipation)]
        r = [solve(testcase(name; eps, T = 0.1, scheme = Scheme(kind = kind))).history.ap_residual[end]
             for eps in (1e-3, 1e-4, 1e-5)]
        @test all(95 .< r[1:2] ./ r[2:3] .< 105)
    end
end

@testset "CSV output" begin
    dir = mktempdir()
    r, paths = run_case(:riemann; dir, N = 400, T = 2.0, snapshots = [1.0])
    @test length(paths) == 3 && all(isfile, paths)
    meta, names, data = read_csv(paths[2])
    @test meta["case"] == "riemann" && meta["N"] == "400" && meta["T"] == "2.0"
    @test haskey(meta, "shock_speed")
    @test names == ["x", "rho", "u", "phi", "rho_exact", "u_exact", "phi_exact"]
    @test data[:, 2] == r.rho
    meta, names, data = read_csv(paths[3])
    @test names[1:3] == ["t", "dt", "energy"] && size(data, 1) == r.steps + 1
end
