# Convergence study for one test case on a sequence of refined grids.
#
#   julia --project=scripts scripts/convergence.jl <case | file.toml> [--scheme shift|dissipation|second] [--lambda λ[,λ2,...]]
#         [--T T]
#         [--N 100,200,400,800] [--dtmax c] [--fixed c] [--var rho,u,phi,m,ap] [--ref Nref]
#         [--limiter minmod|none] [--c_D c] [--eta_c c] [--plot yes|out.png] [--slope 2] [--norm L1,L2,Linf]
#         [--csv yes|out.csv]
#
# The scheme is the file's; --scheme changes its kind (keeping the file's settings when the kind is the same),
# and --limiter, --c_D then set those of a second-order scheme, --eta_c the factor of η in the shift scheme,
# η_e = c (ρ̄ⁿ_e)²/ρⁿʰ_e (default 5/4, (43); file names then get _eta<c>). Unknown options are an error.
# --plot draws the L1, L2 and L∞ errors against N on log–log axes, one panel per variable (no title), with black
# reference lines of the slopes in --slope (default: the scheme's order, 1 or 2; e.g. 1,2), solid for the first
# slope and dashed, dotted for the next ones; style of scripts/style.jl. --norm keeps only those norms in the plot
# (default all three; the tables always show all). `yes` writes convergence_<case>_<scheme>[_<limiter>]_lambda<λ>.png
# in the current directory, with _<norms> appended when --norm leaves out a norm (e.g. _L1).
# Errors are the normalized norms (Δx/|Ω|)Σ|e|, ((Δx/|Ω|)Σe²)^½ and max|e|, so L1 ≤ L2 ≤ L∞.
# One table per variable (default: all of ρ, u, φ and m = ρu); each grid is solved once. The variable ap is the
# quasineutral (AP) residual e^φ - ρ, measured against 0: its norms, which should be O(λ²) on every grid.
# --lambda with several values repeats the study for each λ; the plot then has one curve per λ (and norm) in each
# panel, and `yes` names it with all of them (lambda1.0_0.001_1.0e-6).
# --csv writes the errors, for scripts/plot_convergence.py (the same plot with matplotlib): `# key = value` header
# lines (case, scheme, limiter, c_D or eta_c, T, timestep, reference), then one row per (λ, N) with the columns lambda, N and
# <var>_L1, <var>_L2, <var>_Linf of each variable. `yes` writes convergence_<case>_<scheme>[_<limiter>]_lambda<λ>.csv.
# --dtmax c caps the scheme's Δt rule at c·Δx; --fixed c uses Δt = c·Δx (§5.4 order studies: 0.4).
# Errors are measured against the exact solution when the case has one; otherwise against a
# reference computed by the same scheme on Nref cells (default 4·max N), restricted by cell averaging.
using EPB_AP, Printf
"--plot" in ARGS && (using Plots; include("style.jl"))

const OPTIONS = ("scheme", "lambda", "T", "N", "dtmax", "fixed", "var", "ref", "limiter", "c_D", "eta_c", "plot",
                 "slope", "norm", "csv")
const NORMS = ("L1", "L2", "Linf")

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

# Scheme of the case file, with the kind from --scheme, the second-order settings from --limiter, --c_D and the
# factor of η of the shift scheme from --eta_c.
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
    if haskey(o, "eta_c")
        s isa Scheme && s.kind == :shift || error("--eta_c needs the shift scheme (--scheme shift)")
        s = Scheme(kind = :shift, eta_c = parse(Float64, o["eta_c"]))
    end
    return s
end

field(r, var) = var == "ap" ? exp.(r.phi) .- r.rho : field(r.rho, r.u, r.phi, var)
field(ρ, u, φ, var) = var == "rho" ? ρ : var == "u" ? u : var == "phi" ? φ : ρ .* u
restrict(v, N) = vec(sum(reshape(v, length(v) ÷ N, N); dims = 1)) ./ (length(v) ÷ N)

function main(args)
    name, o = parse_args(args)
    kw = Dict{Symbol,Any}()
    kw[:scheme] = scheme_option(testcase(name).scheme, o)
    haskey(o, "T") && (kw[:T] = parse(Float64, o["T"]))
    haskey(o, "dtmax") && (kw[:timestep] = TimeStep(dtmax = parse(Float64, o["dtmax"])))
    haskey(o, "fixed") && (kw[:timestep] = TimeStep(mode = :fixed, cfl = parse(Float64, o["fixed"])))
    Ns = parse.(Int, split(get(o, "N", "100,200,400,800"), ","))
    vars = split(get(o, "var", "rho,u,phi,m"), ",")
    all(in(("rho", "u", "phi", "m", "ap")), vars) || error("--var takes a list of rho, u, phi, m, ap")
    haskey(o, "slope") && !haskey(o, "plot") && error("--slope needs --plot")
    haskey(o, "norm") && !haskey(o, "plot") && error("--norm needs --plot")
    norms = split(get(o, "norm", join(NORMS, ",")), ",")
    all(in(NORMS), norms) || error("--norm takes a list of $(join(NORMS, ", "))")
    lambdas = haskey(o, "lambda") ? parse.(Float64, split(o["lambda"], ",")) : [testcase(name).lambda]
    tables = [study(name, merge(kw, Dict(:lambda => λ)), Ns, vars, o) for λ in lambdas]
    if haskey(o, "plot")
        case = testcase(name; kw...)
        out = o["plot"] == "yes" ? default_plot_name(case, lambdas, norms) : o["plot"]
        slopes = parse.(Float64, split(get(o, "slope", case.scheme isa SecondOrder ? "2" : "1"), ","))
        plot_errors(Ns, vars, lambdas, tables, slopes, norms, out)
    end
    if haskey(o, "csv")
        case = testcase(name; kw...)
        out = o["csv"] == "yes" ? replace(default_plot_name(case, lambdas, NORMS), r"\.png$" => ".csv") : o["csv"]
        write_csv(out, case, Ns, vars, lambdas, tables, o)
    end
end

# The errors of all (λ, N), one row each, after the settings as `# key = value` lines.
function write_csv(out, case, Ns, vars, lambdas, tables, o)
    s = case.scheme
    open(out, "w") do io
        println(io, "# case = ", case.name)
        println(io, "# scheme = ", EPB_AP.scheme_name(s))
        s isa SecondOrder && println(io, "# limiter = ", s.limiter, "\n# c_D = ", s.c_D)
        s isa Scheme && s.kind == :shift && println(io, "# eta_c = ", s.eta_c)
        println(io, "# T = ", case.T)
        println(io, "# timestep = ", haskey(o, "fixed") ? "fixed $(o["fixed"])" :
                                     haskey(o, "dtmax") ? "cfl, dtmax $(o["dtmax"])" : "cfl")
        println(io, "# reference = ", case.exact === nothing ? "N = $(get(o, "ref", string(4maximum(Ns))))" : "exact")
        println(io, join(["lambda"; "N"; [v * "_" * n for v in vars for n in NORMS]], ","))
        for (λ, table) in zip(lambdas, tables), (i, N) in enumerate(Ns)
            println(io, join([string(λ); string(N); [string(e) for v in vars for e in table[v][i]]], ","))
        end
    end
    println("\nwrote ", out)
end

# Errors of the case with settings kw on the grids Ns: prints one table per variable and returns them,
# variable => [(L1, L2, L∞) per N].
function study(name, kw, Ns, vars, o)
    case = testcase(name; kw...)
    run(N) = (r = solve(testcase(name; kw..., N)); r.status == :completed || error("run with N = $N: $(r.status)"); r)
    if case.exact === nothing
        Nref = parse(Int, get(o, "ref", string(4maximum(Ns))))
        all(Nref .% Ns .== 0) || error("Nref = $Nref must be a multiple of every N")
        rref = run(Nref)
        reference = (r, var) -> var == "ap" ? zero(r.rho) : restrict(field(rref, var), length(r.rho))
        println("Reference: same scheme on N = $Nref cells, restricted by cell averaging (no exact solution).")
    else
        reference = (r, var) -> var == "ap" ? zero(r.rho) : field(case.exact(r.grid.x, r.t)..., var)
        println("Reference: exact solution.")
    end
    s = case.scheme
    header = @sprintf("case = %s, scheme = %s%s, λ = %g, T = %g", case.name, EPB_AP.scheme_name(s),
                      s isa SecondOrder ? " (limiter = $(s.limiter), c_D = $(s.c_D))" :
                      s.kind == :shift ? " (eta_c = $(s.eta_c))" : "", case.lambda, case.T)
    println(header)
    runs = map(run, Ns)
    table = Dict{String,Any}()
    for var in vars
        @printf("\nvariable = %s\n", var)
        @printf("%6s  %11s %6s  %11s %6s  %11s %6s\n", "N", "L1", "order", "L2", "order", "Linf", "order")
        prev = nothing
        for (N, r) in zip(Ns, runs)
            e = field(r, var) .- reference(r, var)
            w = 1 / length(e)                                     # Δx/|Ω|: normalized norms, L1 ≤ L2 ≤ L∞
            errs = (w * sum(abs, e), sqrt(w * sum(abs2, e)), maximum(abs, e))
            ord = prev === nothing ? ("-", "-", "-") :
                  map((a, b) -> @sprintf("%.2f", log(a / b) / log(N / prev[1])), prev[2], errs)
            @printf("%6d  %11.3e %6s  %11.3e %6s  %11.3e %6s\n", N, errs[1], ord[1], errs[2], ord[2], errs[3], ord[3])
            prev = (N, errs)
            push!(get!(table, var, []), errs)
        end
    end
    println()
    return table
end

function default_plot_name(case, lambdas, norms)
    s = case.scheme
    lim = s isa SecondOrder ? "_$(s.limiter)" : ""
    nrm = length(norms) < length(NORMS) ? "_" * join(norms, "_") : ""
    return "convergence_$(case.name)_$(EPB_AP.scheme_name(s))$(lim)$(EPB_AP.eta_tag(s))_lambda$(join(lambdas, "_"))$(nrm).png"
end

# λ in a legend: 1, 0.5, 10⁻³, 2.5·10⁻⁶.
function lambda_label(λ)
    1e-2 <= λ < 1e4 && return "λ = " * @sprintf("%g", λ)
    e = floor(Int, log10(λ)); m = round(λ / 10.0^e; sigdigits = 6)
    sup = Dict(zip("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    p = "10" * map(c -> sup[c], string(e))
    return "λ = " * (m == 1 ? p : "$(m)·$p")
end

# Log–log errors against N in the given norms, one panel per variable and one curve per (λ, norm), with reference
# slopes just below the lowest first point.
function plot_errors(Ns, vars, lambdas, tables, slopes, norms, out)
    panels = map(vars) do var
        p = plot(; xscale = :log10, yscale = :log10, xlabel = "N", ylabel = "error in $(symbol(var))",
                 xticks = (Ns, string.(Ns)), legend = :topright)
        curves = [(λ, name, [errs[findfirst(==(name), NORMS)] for errs in table[var]])
                  for (λ, table) in zip(lambdas, tables) for name in norms]
        if var == "ap"                                    # no order: the level λ² instead
            for λ in lambdas
                reference!(p, Ns, fill(λ^2, length(Ns)); label = length(lambdas) == 1 ? "λ²" : "λ², $(lambda_label(λ))",
                           ls = :dot)
            end
        else
            y0 = 0.5minimum(c[3][1] for c in curves)
            for (j, q) in enumerate(slopes)
                reference!(p, Ns, y0 .* (Ns[1] ./ Ns) .^ q; label = "slope $(q == round(q) ? Int(q) : q)",
                           ls = (:solid, :dash, :dot, :dashdot)[mod1(j, 4)])
            end
        end
        for (k, (λ, name, e)) in enumerate(curves)
            nl = name == "Linf" ? "L∞" : name
            label = length(lambdas) == 1 ? nl : length(norms) == 1 ? lambda_label(λ) : "$nl, $(lambda_label(λ))"
            curve!(p, Ns, e, k; label)
        end
        p
    end
    n = length(panels)
    nc = n == 4 ? 2 : min(n, 3); nr = cld(n, nc)
    fig = plot(panels...; layout = (nr, nc), size = (420nc, 400nr), margin = 4Plots.mm, left_margin = 10Plots.mm, bottom_margin = 8Plots.mm)
    savefig(fig, out)
    println("\nwrote ", out)
end

main(ARGS)
