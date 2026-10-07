"""
    TestCase

A test problem of SPEC §7: periodic domain, initial data `init(x) -> (ρ, u)`, default λ, T, N,
scheme and time-step rule, an optional `exact(x, t) -> (ρ, u, φ)`, snapshot times and
`diagnostics(case, result) -> Vector{Pair}` written to the data file header.
"""
Base.@kwdef struct TestCase
    name::Symbol
    domain::Tuple{Float64,Float64}
    init::Function
    lambda::Float64
    T::Float64
    N::Int
    scheme::Scheme = Scheme()
    timestep::TimeStep = TimeStep()
    exact::Union{Nothing,Function} = nothing
    snapshots::Vector{Float64} = Float64[]
    diagnostics::Function = (case, res) -> Pair{String,Any}[]
end

const CASES = (:smooth, :steep, :nearvacuum, :bump, :step, :expansion, :soliton, :ap, :simplewave, :riemann)

const DISSIPATION = Scheme(kind = :dissipation)

function case_defaults(name::Symbol)
    unit = (0.0, 1.0)
    name == :smooth && return (domain = unit, lambda = 1.0, T = 0.5, N = 200,
        init = x -> (1 .+ 0.3sinpi.(2x), 0.2cospi.(2x)))
    name == :steep && return (domain = unit, lambda = 1.0, T = 0.5, N = 200,
        init = x -> (1 .+ 0.8sinpi.(2x) .^ 3, 0.25cospi.(2x)))
    name == :nearvacuum && return (domain = unit, lambda = 1.0, T = 0.5, N = 200,
        init = x -> (1e-4 .+ (1 - 1e-4) .* exp.(-60 .* (x .- 0.5) .^ 2), zero(x)))
    name == :bump && return (domain = unit, lambda = 1.0, T = 2.0, N = 200, scheme = DISSIPATION,
        init = x -> (1 .+ 0.5exp.(-(x .- 0.5) .^ 2 ./ 0.005), one.(x)))
    name == :step && return (domain = unit, lambda = 0.1, T = 1.0, N = 200, scheme = DISSIPATION,
        init = x -> (ifelse.(0.25 .< x .< 0.75, 1.0, 0.5), zero(x)))
    name == :expansion && return (domain = unit, lambda = 1.0, T = 0.6, N = 200, scheme = DISSIPATION,
        init = x -> (one.(x), sinpi.(2x)))
    if name == :soliton
        c, a, b = 1.3, 0.0, 40.0
        exact = soliton_exact(SolitonProfile(c, 1.0, (b - a) / 2), a, b)
        return (domain = (a, b), lambda = 1.0, T = (b - a) / c, N = 400, scheme = DISSIPATION,
            exact = exact, init = x -> exact(x, 0.0)[1:2])
    end
    name == :ap && return (domain = unit, lambda = 1.0, T = 0.2, N = 200, scheme = DISSIPATION,
        timestep = TimeStep(mode = :fixed, cfl = 0.25),
        init = x -> (1 .+ 0.2sinpi.(2x), 0.2cospi.(2x)))
    name == :simplewave && return (domain = unit, lambda = 1e-4, T = 0.3, N = 200, scheme = DISSIPATION,
        timestep = TimeStep(dtmax = 0.4), exact = simplewave_exact,
        init = x -> simplewave_exact(x, 0.0)[1:2])
    name == :riemann && return (domain = (-80.0, 100.0), lambda = 1e-4, T = 8.0, N = 2000,
        snapshots = [4.0], exact = riemann_exact, init = x -> riemann_exact(x, 0.0)[1:2],
        diagnostics = riemann_diagnostics)
    throw(ArgumentError("unknown test case $name; available: $(join(CASES, ", "))"))
end

"""
    testcase(name; kwargs...)

Test case `name` (one of `EPB_AP.CASES`) with defaults from SPEC §7; any field of
[`TestCase`](@ref) (e.g. `lambda`, `T`, `N`, `scheme`, `timestep`) can be overridden.
"""
testcase(name::Symbol; kw...) = TestCase(; name, case_defaults(name)..., kw...)

"Run a test case: `solve` on its grid with its defaults."
function solve(case::TestCase)
    g = Grid(case.domain..., case.N)
    ρ0, u0 = case.init(g.x)
    return solve(ρ0, u0, g, case.lambda, case.T; scheme = case.scheme, timestep = case.timestep,
                 snapshots = case.snapshots)
end

"Shock speed (x_s(T) - x_s(T₁))/(T - T₁) and its relative error (§5.5)."
function riemann_diagnostics(case, res)
    isempty(res.snapshots) && return Pair{String,Any}[]
    t1, ρ1 = res.snapshots[1][1], res.snapshots[1][2]
    g = res.grid
    s = (shock_position(g.x, res.rho, res.t, g.dx) - shock_position(g.x, ρ1, t1, g.dx)) / (res.t - t1)
    us = riemann_speeds(0.5)[2]
    return Pair{String,Any}["shock_speed" => s, "shock_speed_exact" => us,
                            "shock_speed_rel_error" => abs(s - us) / us]
end
