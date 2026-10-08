using EPB_AP: CASES, lap, riemann_speeds, SolitonProfile

@testset "test case definitions" begin
    for name in CASES
        c = testcase(name)
        @test c.name == (name == :riemann_speed ? :riemann : name)
        g = Grid(c.domain..., 64)
        ρ, u = c.init(g.x)
        @test length(ρ) == length(u) == 64 && all(ρ .> 0)
    end
    @test testcase(:smooth; lambda = 1e-2, N = 50).lambda == 1e-2
    @test_throws ArgumentError testcase(:nope)
    # parameter files: the values of SPEC §7, and a user file of the same form
    c = testcase(:soliton)
    @test c.domain == (0.0, 40.0) && c.N == 400 && c.T ≈ 40 / 1.3 && c.scheme == Scheme(kind = :dissipation)
    @test testcase(:ap).timestep == TimeStep(mode = :fixed, cfl = 0.25)
    file = joinpath(mktempdir(), "mine.toml")
    write(file, replace(read(joinpath(EPB_AP.PARAMS, "step.toml"), String), "N = 200" => "N = 64"))
    c = testcase(file)
    @test c.name == :step && c.N == 64 && c.lambda == 0.1 && c.scheme.kind == :dissipation
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

@testset "asymptotic preservation: residual ∝ λ²" begin
    for (name, kind) in [(:smooth, :shift), (:ap, :dissipation)]
        r = [solve(testcase(name; lambda, T = 0.1, scheme = Scheme(kind = kind))).history.ap_residual[end]
             for lambda in (1e-3, 1e-4, 1e-5)]
        @test all(95 .< r[1:2] ./ r[2:3] .< 105)
    end
end

@testset "Riemann problem with walls [1, §7.4]: shock position and middle state" begin
    c = testcase(:riemann)
    @test c.bc == :wall && c.N == 400 && c.T == 50
    r = solve(c)
    h = r.history
    @test r.status == :completed && abs(r.u[1]) < 1e-3 && abs(r.u[end]) < 1e-3   # only tails reach the walls
    @test maximum(abs, h.mass .- h.mass[1]) <= 1e-13 * h.mass[1]
    @test all(diff(h.energy) .<= 1e3 * eps() * max(1, abs(h.energy[1])))
    d = Dict(c.diagnostics(c, r))
    @test d["shock_position_exact"] ≈ 58.89 atol = 0.01               # [1]: 58.9
    @test abs(d["shock_position"] - d["shock_position_exact"]) < 0.5   # Δx = 0.45
    @test r.u[argmin(abs.(r.grid.x))] ≈ riemann_speeds(0.5)[1] atol = 1e-3
end

@testset "CSV output" begin
    dir = mktempdir()
    paths = run_case(:riemann; dir, N = 400, T = 2.0, snapshots = [1.0])
    @test paths isa Vector{String}
    r = solve(testcase(:riemann; N = 400, T = 2.0, snapshots = [1.0]))
    @test length(paths) == 3 && all(isfile, paths)
    meta, names, data = read_csv(paths[2])
    @test meta["case"] == "riemann" && meta["N"] == "400" && meta["T"] == "2.0" && meta["bc"] == "wall"
    @test haskey(meta, "shock_speed")
    @test names == ["x", "rho", "u", "phi", "rho_init", "u_init", "phi_init", "rho_exact", "u_exact", "phi_exact"]
    @test data[:, 2] == r.rho
    @test data[:, 5:6] == hcat(testcase(:riemann).init(r.grid.x)...)
    meta, names, data = read_csv(paths[3])
    @test names[1:3] == ["t", "dt", "energy"] && size(data, 1) == r.steps + 1
end

@testset "second order on every case: mass, positivity, AP, order" begin
    for name in EPB_AP.CASES
        c = testcase(name; scheme = SecondOrder())
        r = solve(testcase(name; scheme = SecondOrder(), T = min(c.T, 0.1 * (c.domain[2] - c.domain[1]))))
        h = r.history
        @test r.status == :completed
        @test maximum(abs, h.mass .- h.mass[1]) <= 1e-13 * h.mass[1]
        @test all(h.min_rho .> 0)
    end
    r = [solve(testcase(:ap; lambda, T = 0.1, scheme = SecondOrder())).history.ap_residual[end]
         for lambda in (1e-3, 1e-4, 1e-5)]
    @test all(95 .< r[1:2] ./ r[2:3] .< 105)
    errs = map((100, 200, 400)) do N
        c = testcase(:simplewave; N, scheme = SecondOrder(), timestep = TimeStep(mode = :fixed, cfl = 0.4))
        r = solve(c)
        r.grid.dx * sum(abs, r.rho .- c.exact(r.grid.x, r.t)[1])
    end
    @test all(1.9 .< log2.(errs[1:2] ./ errs[2:3]) .< 2.1)
end

@testset "second order from a parameter file and in output names" begin
    dir = mktempdir()
    file = joinpath(dir, "step2.toml")
    write(file, replace(read(joinpath(EPB_AP.PARAMS, "step.toml"), String),
                        r"\[scheme\][^\[]*" => "[scheme]\nkind = \"second\"\nlimiter = \"none\"\nc_D = 2.0\n\n"))
    c = testcase(file; T = 0.05)
    @test c.scheme == SecondOrder(limiter = :none, c_D = 2.0)
    paths = run_case(file; dir, T = 0.05)
    @test occursin("step_second_N200", paths[1])
    meta, _, _ = read_csv(paths[1])
    @test meta["scheme"] == "second" && meta["limiter"] == "none"
    @test_throws ArgumentError EPB_AP.scheme_from(Dict("kind" => "third"))
end
