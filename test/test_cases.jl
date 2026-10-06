using EPB_AP: CASES

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

@testset "properties along runs: mass, positivity, energy" begin
    for kind in (:shift, :dissipation), ε in (1.0, 1e-4)
        r = solve(testcase(:smooth; eps = ε, T = 0.2, scheme = Scheme(kind = kind)))
        h = r.history
        @test r.status == :completed && r.t == 0.2
        @test maximum(abs, h.mass .- h.mass[1]) <= 1e-13 * h.mass[1]
        @test all(h.min_rho .> 0)
        @test all(diff(h.energy) .<= 1e3 * eps() * max(1, abs(h.energy[1])))
        @test all(filter(!isnan, h.cfl_ratio) .<= 1)
    end
end

@testset "asymptotic preservation: residual ∝ ε²" begin
    r = [solve(testcase(:smooth; eps, T = 0.1)).history.ap_residual[end] for eps in (1e-3, 1e-4, 1e-5)]
    @test all(95 .< r[1:2] ./ r[2:3] .< 105)
end

@testset "CSV output" begin
    dir = mktempdir()
    r, paths = run_case(:smooth; dir, N = 50, T = 0.1, snapshots = [0.05])
    @test length(paths) == 3 && all(isfile, paths)
    meta, names, data = read_csv(paths[2])
    @test meta["case"] == "smooth" && meta["N"] == "50" && meta["T"] == "0.1" && meta["eps"] == "1.0"
    @test names == ["x", "rho", "u", "phi"]
    @test data[:, 2] == r.rho
    @test read_csv(paths[1])[1]["t"] == "0.05"
    meta, names, data = read_csv(paths[3])
    @test names[1:3] == ["t", "dt", "energy"] && size(data, 1) == r.steps + 1
end
