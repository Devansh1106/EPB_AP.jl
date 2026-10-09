# Convergence study for one test case on a sequence of refined grids.
#
#   julia --project=scripts scripts/convergence.jl <case | file.toml> [--scheme shift|dissipation|second] [--lambda λ] [--T T]
#         [--N 100,200,400,800] [--dtmax c] [--fixed c] [--var rho,u,phi,m] [--ref Nref]
#         [--limiter minmod|none] [--c_D c]
#
# The scheme is the file's; --scheme changes its kind (keeping the file's settings when the kind is the same),
# and --limiter, --c_D then set those of a second-order scheme. Unknown options are an error.
# One table per variable (default: all of ρ, u, φ and m = ρu); each grid is solved once.
# --dtmax c caps the scheme's Δt rule at c·Δx; --fixed c uses Δt = c·Δx (§5.4 order studies: 0.4).
# Errors are measured against the exact solution when the case has one; otherwise against a
# reference computed by the same scheme on Nref cells (default 4·max N), restricted by cell averaging.
using EPB_AP, Printf

const OPTIONS = ("scheme", "lambda", "T", "N", "dtmax", "fixed", "var", "ref", "limiter", "c_D")

function parse_args(args)
    usage = "usage: convergence.jl <case | file.toml> " * join(("[--$k ..]" for k in OPTIONS), " ")
    isempty(args) && error(usage)
    iseven(length(args) - 1) || error("every option needs a value\n" * usage)
    opts = Dict{String,String}()
    for i in 2:2:length(args)
        k = startswith(args[i], "--") ? args[i][3:end] : ""
        k in OPTIONS || error("unknown option $(args[i])\n" * usage)
        opts[k] = args[i+1]
    end
    return endswith(args[1], ".toml") ? args[1] : Symbol(args[1]), opts
end

# Scheme of the case file, with the kind from --scheme and the second-order settings from --limiter, --c_D.
function scheme_option(file_scheme, o)
    s = file_scheme
    if haskey(o, "scheme") && String(EPB_AP.scheme_name(s)) != o["scheme"]
        s = EPB_AP.scheme_from(Dict("kind" => o["scheme"]))
    end
    if haskey(o, "limiter") || haskey(o, "c_D")
        s isa SecondOrder || error("--limiter and --c_D need a second-order scheme (--scheme second)")
        s = SecondOrder(limiter = Symbol(get(o, "limiter", String(s.limiter))),
                        c_D = parse(Float64, get(o, "c_D", string(s.c_D))), cfl = s.cfl)
    end
    return s
end

field(r, var) = var == "rho" ? r.rho : var == "u" ? r.u : var == "phi" ? r.phi : r.rho .* r.u
field(ρ, u, φ, var) = var == "rho" ? ρ : var == "u" ? u : var == "phi" ? φ : ρ .* u
restrict(v, N) = vec(sum(reshape(v, length(v) ÷ N, N); dims = 1)) ./ (length(v) ÷ N)

function main(args)
    name, o = parse_args(args)
    kw = Dict{Symbol,Any}()
    kw[:scheme] = scheme_option(testcase(name).scheme, o)
    haskey(o, "lambda") && (kw[:lambda] = parse(Float64, o["lambda"]))
    haskey(o, "T") && (kw[:T] = parse(Float64, o["T"]))
    haskey(o, "dtmax") && (kw[:timestep] = TimeStep(dtmax = parse(Float64, o["dtmax"])))
    haskey(o, "fixed") && (kw[:timestep] = TimeStep(mode = :fixed, cfl = parse(Float64, o["fixed"])))
    Ns = parse.(Int, split(get(o, "N", "100,200,400,800"), ","))
    vars = split(get(o, "var", "rho,u,phi,m"), ",")
    case = testcase(name; kw...)
    run(N) = (r = solve(testcase(name; kw..., N)); r.status == :completed || error("run with N = $N: $(r.status)"); r)
    if case.exact === nothing
        Nref = parse(Int, get(o, "ref", string(4maximum(Ns))))
        all(Nref .% Ns .== 0) || error("Nref = $Nref must be a multiple of every N")
        rref = run(Nref)
        reference = (r, var) -> restrict(field(rref, var), length(r.rho))
        println("Reference: same scheme on N = $Nref cells, restricted by cell averaging (no exact solution).")
    else
        reference = (r, var) -> field(case.exact(r.grid.x, r.t)..., var)
        println("Reference: exact solution.")
    end
    s = case.scheme
    @printf("case = %s, scheme = %s%s, λ = %g, T = %g\n", case.name, EPB_AP.scheme_name(s),
            s isa SecondOrder ? " (limiter = $(s.limiter), c_D = $(s.c_D))" : "", case.lambda, case.T)
    runs = map(run, Ns)
    for var in vars
        @printf("\nvariable = %s\n", var)
        @printf("%6s  %11s %6s  %11s %6s  %11s %6s\n", "N", "L1", "order", "L2", "order", "Linf", "order")
        prev = nothing
        for (N, r) in zip(Ns, runs)
            e = field(r, var) .- reference(r, var)
            errs = (r.grid.dx * sum(abs, e), sqrt(r.grid.dx * sum(abs2, e)), maximum(abs, e))
            ord = prev === nothing ? ("-", "-", "-") :
                  map((a, b) -> @sprintf("%.2f", log(a / b) / log(N / prev[1])), prev[2], errs)
            @printf("%6d  %11.3e %6s  %11.3e %6s  %11.3e %6s\n", N, errs[1], ord[1], errs[2], ord[2], errs[3], ord[3])
            prev = (N, errs)
        end
    end
end

main(ARGS)
