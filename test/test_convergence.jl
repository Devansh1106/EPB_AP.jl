@testset "observed order of accuracy ≈ 1 (simple wave, exact solution)" begin
    for kind in (:shift, :dissipation)
        errs = map((100, 200, 400)) do N
            c = testcase(:simplewave; N, scheme = Scheme(kind = kind))
            r = solve(c)
            r.grid.dx * sum(abs.(r.rho .- c.exact(r.grid.x, r.t)[1]))
        end
        p = log2.(errs[1:2] ./ errs[2:3])
        @test all(0.9 .< p .< 1.1)
    end
end
