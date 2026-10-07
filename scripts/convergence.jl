# Convergence study for one test case on a sequence of refined grids.
#
#   julia --project=scripts scripts/convergence.jl <case | file.toml> [--scheme shift|dissipation] [--lambda λ] [--T T]
#         [--N 100,200,400,800] [--dtmax c] [--var rho|u|phi|m] [--ref Nref]
#
# Errors are measured against the exact solution when the case has one; otherwise against a
# reference computed by the same scheme on Nref cells (default 4·max N), restricted by cell averaging.
using EPB_AP, Printf

function parse_args(args)
    isempty(args) && error("usage: convergence.jl <case> [--scheme ..] [--lambda ..] [--T ..] [--N ..] [--dtmax ..] [--var ..] [--ref ..]")
    opts = Dict(args[i][3:end] => args[i+1] for i in 2:2:length(args)-1)
    return endswith(args[1], ".toml") ? args[1] : Symbol(args[1]), opts
end

field(r, var) = var == "rho" ? r.rho : var == "u" ? r.u : var == "phi" ? r.phi : r.rho .* r.u
field(ρ, u, φ, var) = var == "rho" ? ρ : var == "u" ? u : var == "phi" ? φ : ρ .* u
restrict(v, N) = vec(sum(reshape(v, length(v) ÷ N, N); dims = 1)) ./ (length(v) ÷ N)

function main(args)
    name, o = parse_args(args)
    kw = Dict{Symbol,Any}()
    haskey(o, "scheme") && (kw[:scheme] = Scheme(kind = Symbol(o["scheme"])))
    haskey(o, "lambda") && (kw[:lambda] = parse(Float64, o["lambda"]))
    haskey(o, "T") && (kw[:T] = parse(Float64, o["T"]))
    haskey(o, "dtmax") && (kw[:timestep] = TimeStep(dtmax = parse(Float64, o["dtmax"])))
    Ns = parse.(Int, split(get(o, "N", "100,200,400,800"), ","))
    var = get(o, "var", "rho")
    case = testcase(name; kw...)
    run(N) = (r = solve(testcase(name; kw..., N)); r.status == :completed || error("run with N = $N: $(r.status)"); r)
    if case.exact === nothing
        Nref = parse(Int, get(o, "ref", string(4maximum(Ns))))
        all(Nref .% Ns .== 0) || error("Nref = $Nref must be a multiple of every N")
        ref = field(run(Nref), var)
        reference = r -> restrict(ref, length(r.rho))
        println("Reference: same scheme on N = $Nref cells, restricted by cell averaging (no exact solution).")
    else
        reference = r -> field(case.exact(r.grid.x, r.t)..., var)
        println("Reference: exact solution.")
    end
    @printf("case = %s, scheme = %s, λ = %g, T = %g, variable = %s\n",
            case.name, case.scheme.kind, case.lambda, case.T, var)
    @printf("%6s  %11s %6s  %11s %6s  %11s %6s\n", "N", "L1", "order", "L2", "order", "Linf", "order")
    prev = nothing
    for N in Ns
        r = run(N)
        e = field(r, var) .- reference(r)
        errs = (r.grid.dx * sum(abs, e), sqrt(r.grid.dx * sum(abs2, e)), maximum(abs, e))
        ord = prev === nothing ? ("-", "-", "-") :
              map((a, b) -> @sprintf("%.2f", log(a / b) / log(N / prev[1])), prev[2], errs)
        @printf("%6d  %11.3e %6s  %11.3e %6s  %11.3e %6s\n", N, errs[1], ord[1], errs[2], ord[2], errs[3], ord[3])
        prev = (N, errs)
    end
end

main(ARGS)
