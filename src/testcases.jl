"""
    TestCase

A test problem of SPEC §7: periodic domain, initial data `init(x) -> (ρ, u)`, default ε, T, N,
scheme and time-step rule, an optional `exact(x, t) -> (ρ, u, φ)`, snapshot times and
`diagnostics(case, result) -> Vector{Pair}` written to the data file header.
"""
Base.@kwdef struct TestCase
    name::Symbol
    domain::Tuple{Float64,Float64}
    init::Function
    eps::Float64
    T::Float64
    N::Int
    scheme::Scheme = Scheme()
    timestep::TimeStep = TimeStep()
    exact::Union{Nothing,Function} = nothing
    snapshots::Vector{Float64} = Float64[]
    diagnostics::Function = (case, res) -> Pair{String,Any}[]
end

const CASES = (:smooth,)

function case_defaults(name::Symbol)
    unit = (0.0, 1.0)
    name == :smooth && return (domain = unit, eps = 1.0, T = 0.5, N = 200,
        init = x -> (1 .+ 0.3sinpi.(2x), 0.2cospi.(2x)))
    throw(ArgumentError("unknown test case $name; available: $(join(CASES, ", "))"))
end

"""
    testcase(name; kwargs...)

Test case `name` (one of `EPB_AP.CASES`) with defaults from SPEC §7; any field of
[`TestCase`](@ref) (e.g. `eps`, `T`, `N`, `scheme`, `timestep`) can be overridden.
"""
testcase(name::Symbol; kw...) = TestCase(; name, case_defaults(name)..., kw...)

"Run a test case: `solve` on its grid with its defaults."
function solve(case::TestCase)
    g = Grid(case.domain..., case.N)
    ρ0, u0 = case.init(g.x)
    return solve(ρ0, u0, g, case.eps, case.T; scheme = case.scheme, timestep = case.timestep,
                 snapshots = case.snapshots)
end
