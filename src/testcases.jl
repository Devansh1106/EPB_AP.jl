"""
    TestCase

A test problem of SPEC §7, built from a parameter file: domain and boundary (`bc`), initial data `init(x) -> (ρ, u)`, default λ, T, N,
scheme and time-step rule, an optional `exact(x, t) -> (ρ, u, φ)`, snapshot times and
`diagnostics(case, result) -> Vector{Pair}` written to the data file header.
"""
Base.@kwdef struct TestCase
    name::Symbol
    domain::Tuple{Float64,Float64}
    bc::Symbol = :periodic
    init::Function
    lambda::Float64
    T::Float64
    N::Int
    scheme::Union{Scheme,SecondOrder} = Scheme()
    timestep::TimeStep = TimeStep()
    exact::Union{Nothing,Function} = nothing
    snapshots::Vector{Float64} = Float64[]
    diagnostics::Function = (case, res) -> Pair{String,Any}[]
end

const CASES = (:smooth, :steep, :nearvacuum, :bump, :step, :expansion, :soliton, :ap, :simplewave,
               :riemann, :riemann_speed)

"Directory of the parameter files `<case>.toml`, one per entry of `CASES`."
const PARAMS = joinpath(dirname(@__DIR__), "params")

# Initial data, exact solution and diagnostics of a case; `p` holds the `[init]` table of its file.
function problem(name::Symbol, domain, λ, p)
    name == :smooth && return (init = x -> (1 .- 0.3sech.(2x), fill(float(p["U"]), length(x))),)
    name == :steep && return (init = x -> (1 .+ 0.8sinpi.(2x) .^ 3, 0.25cospi.(2x)),)
    name == :nearvacuum && return (init = x -> (1e-4 .+ (1 - 1e-4) .* exp.(-60 .* (x .- 0.5) .^ 2), zero(x)),)
    name == :bump && return (init = x -> (1 .+ 0.5exp.(-(x .- 0.5) .^ 2 ./ 0.005), one.(x)),)
    name == :step && return (init = x -> (ifelse.(0.25 .< x .< 0.75, 1.0, 0.5), zero(x)),)
    name == :expansion && return (init = x -> (one.(x), sinpi.(2x)),)
    if name == :soliton
        a, b = domain
        exact = soliton_exact(SolitonProfile(p["c"], λ, (b - a) / 2), a, b)
        return (exact = exact, init = x -> exact(x, 0.0)[1:2])
    end
    name == :ap && return (init = x -> (1 .+ 0.2sinpi.(2x), 0.2cospi.(2x)),)
    name == :simplewave && return (exact = simplewave_exact, init = x -> simplewave_exact(x, 0.0)[1:2])
    if name == :riemann
        exact = (x, t) -> riemann_exact(x, t; nr = p["nr"])
        return (exact = exact, init = x -> exact(x, 0.0)[1:2], diagnostics = riemann_diagnostics(p["nr"]))
    end
    throw(ArgumentError("unknown test case $name; available: $(join(CASES, ", "))"))
end

symbols(d) = (; (Symbol(k) => (v isa String ? Symbol(v) : v) for (k, v) in d)...)

"""
    scheme_from(d)

Scheme from the `[scheme]` table of a parameter file (or any `Dict`): `kind = "second"` gives
[`SecondOrder`](@ref), `"shift"` or `"dissipation"` a first-order [`Scheme`](@ref); the other keys are passed on.
"""
function scheme_from(d)
    get(d, "kind", "shift") == "second" || return Scheme(; symbols(d)...)
    return SecondOrder(; symbols(filter(p -> p.first != "kind", d))...)
end

"""
    testcase(name::Symbol; kwargs...)
    testcase(file::AbstractString; kwargs...)

Test case read from its parameter file, `params/<name>.toml` or any `file` of the same form
(SPEC §7). Keywords override fields of [`TestCase`](@ref), e.g. `lambda`, `T`, `N`, `scheme`.
"""
function testcase(name::Symbol; kw...)
    name in CASES || throw(ArgumentError("unknown test case $name; available: $(join(CASES, ", "))"))
    return testcase(joinpath(PARAMS, "$name.toml"); kw...)
end

function testcase(file::AbstractString; kw...)
    p = TOML.parsefile(file)
    name, domain = Symbol(p["case"]), Tuple(float.(p["domain"]))
    λ = get(kw, :lambda, p["lambda"])
    return TestCase(; name, domain, bc = Symbol(get(p, "bc", "periodic")), lambda = λ, T = p["T"], N = p["N"],
                    snapshots = float.(p["snapshots"]), scheme = scheme_from(p["scheme"]),
                    timestep = TimeStep(; symbols(p["timestep"])...),
                    problem(name, domain, λ, get(p, "init", Dict()))..., kw...)
end

"Run a test case: `solve` on its grid with its defaults."
function solve(case::TestCase)
    g = Grid(case.domain..., case.N; bc = case.bc)
    ρ0, u0 = case.init(g.x)
    return solve(ρ0, u0, g, case.lambda, case.T; scheme = case.scheme, timestep = case.timestep,
                 snapshots = case.snapshots)
end

"""
Shock position x_s(T) against u_s T ([1, Fig. 8]) and, with a snapshot at T₁, the shock speed
(x_s(T) - x_s(T₁))/(T - T₁) and its relative error (§5.5).
"""
riemann_diagnostics(nr) = function (case, res)
    g, us = res.grid, riemann_speeds(nr)[2]
    xs = shock_position(g.x, res.rho, res.t, g.dx; nr)
    d = Pair{String,Any}["shock_position" => xs, "shock_position_exact" => us * res.t]
    isempty(res.snapshots) && return d
    t1, ρ1 = res.snapshots[1][1], res.snapshots[1][2]
    s = (xs - shock_position(g.x, ρ1, t1, g.dx; nr)) / (res.t - t1)
    return [d; Pair{String,Any}["shock_speed" => s, "shock_speed_exact" => us,
                                "shock_speed_rel_error" => abs(s - us) / us]]
end
