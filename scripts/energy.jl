# Energy stability study: the discrete energy (14) along runs of one test case, for several schemes and λ.
#
#   julia --project=scripts scripts/energy.jl <case | file.toml> [--schemes shift,dissipation,second-none,second-minmod]
#         [--lambda 1,1e-2,1e-4,1e-6] [--N N] [--T T] [--plot yes|out.png] [--relative]
#   julia --project=scripts scripts/energy.jl <run>_history.csv [more_history.csv ...] [--plot out.png] [--relative]
#
# A scheme is shift, dissipation or second, optionally with a limiter (second-none, second-minmod); a kind equal to
# the file's keeps the file's settings. For each scheme and λ the table gives E(T) - E(0), the largest step change
# max_n (E^{n+1} - E^n) and the number of steps with E^{n+1} > E^n beyond round-off (1e3·eps·max(1, |E^0|), as in
# the tests). Energy stability is max ΔE ≤ 0, i.e. no increasing step. First order: Theorems 3.5, 3.16; second
# order: no fully discrete inequality is claimed (§4.13), so the table is an observation.
# --plot draws E(t) - E(0) against t, one panel per scheme (FOPD, FODD, SODD; named in its legend; no title), one
# curve per λ, in the style of scripts/style.jl. `yes` writes energy_<case>_N<N>.png in the current directory.
# --relative uses E(t)/E(0) instead of E(t) - E(0), in the table and the plot (the energy (14) itself, unshifted: where
# E(0) < 0, a decreasing energy gives a ratio above 1).
#
# Given history files (written by run_case) instead of a case, nothing is run: the table is that of these runs, and
# one plot shows them all (legend FOPD, FODD, SODD), with the energy of the exact solution, if the case has one, on the
# grid of the first file, as the black reference. Default output: energy_<case>_N<N>_compare.png
# (energy_relative_... with --relative).
using EPB_AP, Printf
("--plot" in ARGS || any(endswith(".csv"), ARGS)) && (using Plots; include("style.jl"))

const OPTIONS = ("schemes", "lambda", "N", "T", "plot")
const CASE_OPTIONS = ("schemes", "lambda", "N", "T")

# Returns (case name or parameter file, or history files), options, relative.
function parse_args(args)
    usage = "usage: energy.jl <case | file.toml | history.csv ...> " * join(("[--$k ..]" for k in OPTIONS), " ") *
            " [--relative]"
    rel = "--relative" in args
    args = filter(!=("--relative"), args)
    n = something(findfirst(startswith("--"), args), length(args) + 1) - 1
    n >= 1 || error(usage)
    srcs, rest = args[1:n], args[n+1:end]
    iseven(length(rest)) || error("every option needs a value\n" * usage)
    opts = Dict{String,String}()
    for i in 1:2:length(rest)
        k = startswith(rest[i], "--") ? rest[i][3:end] : ""
        k in OPTIONS || error("unknown option $(rest[i])\n" * usage)
        opts[k] = rest[i+1]
    end
    if all(endswith(".csv"), srcs)
        bad = intersect(keys(opts), CASE_OPTIONS)
        isempty(bad) || error("--$(join(bad, ", --")) set(s) up runs; history files are already run\n" * usage)
        return srcs, opts, rel
    end
    n == 1 || error("give one case or parameter file, or history files (.csv)\n" * usage)
    return endswith(srcs[1], ".toml") ? srcs[1] : Symbol(srcs[1]), opts, rel
end

measure(E, rel) = rel ? E ./ E[1] : E .- E[1]
measure_name(rel) = rel ? "E(t) / E(0)" : "E(t) − E(0)"

# Scheme `kind[-limiter]`, keeping the file's settings when the kind is the file's.
function scheme_named(file_scheme, name)
    kind, lim = occursin("-", name) ? String.(split(name, "-"; limit = 2)) : (String(name), nothing)
    kind in ("shift", "dissipation", "second") || error("unknown scheme $name")
    s = String(EPB_AP.scheme_name(file_scheme)) == kind ? file_scheme : EPB_AP.scheme_from(Dict("kind" => kind))
    if lim !== nothing
        s isa SecondOrder || error("a limiter needs a second-order scheme: $name")
        s = SecondOrder(limiter = Symbol(lim), c_D = s.c_D, cfl = s.cfl)
    end
    return s
end

table_header(first, rel) = @printf("%-14s %8s  %12s  %12s  %7s  %s\n", first, "λ", rel ? "E(T)/E(0)" : "E(T)-E(0)",
                                   "max ΔE", "#ΔE>0", "status")

function table_row(n, λ, E, status, rel)
    dE = diff(E)
    tol = 1e3 * eps() * max(1, abs(E[1]))
    @printf("%-14s %8g  %12.4e  %12.4e  %7d  %s\n", n, λ, measure(E, rel)[end], isempty(dE) ? NaN : maximum(dE),
            count(>(tol), dE), status)
end

function main(args)
    src, o, rel = parse_args(args)
    src isa Vector && return compare_files(src, o, rel)
    name = src
    kw = Dict{Symbol,Any}()
    haskey(o, "N") && (kw[:N] = parse(Int, o["N"]))
    haskey(o, "T") && (kw[:T] = parse(Float64, o["T"]))
    case = testcase(name; kw...)
    names = split(get(o, "schemes", "shift,dissipation,second-none,second-minmod"), ",")
    schemes = [n => scheme_named(case.scheme, n) for n in names]
    λs = parse.(Float64, split(get(o, "lambda", "1,1e-2,1e-4,1e-6"), ","))
    header = @sprintf("case = %s, N = %d, T = %g", case.name, case.N, case.T)
    println(header, "\n")
    table_header("scheme", rel)
    curves = Dict{String,Any}()
    for (n, s) in schemes, λ in λs
        r = solve(testcase(name; kw..., lambda = λ, scheme = s))
        E = r.history.energy
        table_row(n, λ, E, r.status == :completed ? "completed" : "$(r.status) at t = $(r.t)", rel)
        push!(get!(curves, n, []), (λ, r.history.t, measure(E, rel)))
    end
    if haskey(o, "plot")
        out = o["plot"] == "yes" ? "energy_$(rel ? "relative_" : "")$(case.name)_N$(case.N).png" : o["plot"]
        plot_energy(names, curves, out, rel)
    end
end

# Legend name of a scheme option: FOPD, FODD, SODD (scripts/style.jl), the limiter in brackets.
function legend_name(n)
    kind, lim = occursin("-", n) ? split(n, "-"; limit = 2) : (n, nothing)
    return scheme_label(kind) * (lim === nothing ? "" : " ($lim)")
end

# Legend position away from the curves, which leave (0, E(0)) going down (ΔE ≤ 0) or, for a ratio with E(0) < 0, up.
legend_side(ys) = any(y -> y[end] > y[1], ys) ? :topleft : :bottomleft

# E(t) - E(0) (or E(t)/E(0)) against t, one panel per scheme, one curve per λ.
function plot_energy(names, curves, out, rel)
    panels = map(names) do n
        p = plot(; xlabel = "t", ylabel = measure_name(rel), legend = legend_side([y for (_, _, y) in curves[n]]),
                 legend_title = legend_name(n), legend_title_font_pointsize = 8)
        for (k, (λ, t, y)) in enumerate(curves[n])
            curve!(p, t, y, k; label = "λ = $(λ)")
        end
        p
    end
    n = length(panels)
    nc = n == 4 ? 2 : min(n, 3); nr = cld(n, nc)
    fig = plot(panels...; layout = (nr, nc), size = (420nc, 380nr), margin = 4Plots.mm, left_margin = 10Plots.mm, bottom_margin = 8Plots.mm)
    savefig(fig, out)
    println("\nwrote ", out)
end

# Legend names of runs from their headers: the scheme (FOPD, FODD, SODD), with the limiter, N, λ, T added in turn
# while two runs still share a name; the file name if that is not enough.
function run_labels(paths, metas)
    labs = [scheme_label(get(m, "scheme", "?")) for m in metas]
    for k in ("limiter", "N", "lambda", "T")
        allunique(labs) && break
        length(unique(get(m, k, "") for m in metas)) > 1 || continue
        labs = [l * ", $(k == "lambda" ? "λ" : k)=$(get(m, k, "?"))" for (l, m) in zip(labs, metas)]
    end
    return allunique(labs) ? labs : basename.(paths)
end

# Energy (14) of the exact solution of the run's case, sampled on its grid at (at most 400 of) the times `t`; nothing
# if the case has no exact solution. The case's [init] data are those of params/<case>.toml.
function exact_energy(m, t)
    name = Symbol(get(m, "case", ""))
    name in EPB_AP.CASES || return nothing
    domain = Tuple(parse.(Float64, split(strip(m["domain"], ['[', ']']), ",")))
    case = testcase(name; N = parse(Int, m["N"]), lambda = parse(Float64, m["lambda"]), domain, bc = Symbol(m["bc"]))
    case.exact === nothing && return nothing
    g = Grid(case.domain..., case.N; bc = case.bc)
    ts = length(t) <= 400 ? t : collect(range(first(t), last(t); length = 400))
    E = [energy(case.exact(g.x, s)..., case.lambda, g.dx; wall = EPB_AP.neumann_phi(case.bc)) for s in ts]
    return ts, E
end

function compare_files(paths, o, rel)
    runs = map(paths) do p
        m, names, data = read_csv(p)
        j, it = findfirst(==("energy"), names), findfirst(==("t"), names)
        (j === nothing || it === nothing) && error("$p has no t and energy columns: give a run's _history.csv file")
        (m, data[:, it], data[:, j])
    end
    metas = first.(runs)
    labels = run_labels(paths, metas)
    m1 = metas[1]
    println("case = $(get(m1, "case", "?")), N = $(get(m1, "N", "?")), T = $(get(m1, "T", "?"))\n")
    table_header("run", rel)
    for ((m, t, E), lab) in zip(runs, labels)
        table_row(lab, parse(Float64, get(m, "lambda", "NaN")), E,
                  get(m, "status", "?") == "completed" ? "completed" : "$(get(m, "status", "?")) at t = $(t[end])", rel)
    end
    ex = exact_energy(m1, runs[1][2])
    if ex !== nothing && abs(ex[2][1] - runs[1][3][1]) > 1e-2 * abs(runs[1][3][1])
        @warn "the exact solution's energy at t = 0 ($(ex[2][1])) is not the run's ($(runs[1][3][1])): its initial " *
              "data differ from params/$(m1["case"]).toml; the exact energy is left out"
        ex = nothing
    end
    ys = [measure(E, rel) for (_, _, E) in runs]
    p = plot(; xlabel = "t", ylabel = measure_name(rel), legend = legend_side(ys))
    for (k, ((_, t, _), y, lab)) in enumerate(zip(runs, ys, labels))
        curve!(p, t, y, k; label = lab)
    end
    ex === nothing || reference!(p, ex[1], measure(ex[2], rel); label = "Exact")
    out = get(o, "plot", "yes")
    out == "yes" && (out = "energy_$(rel ? "relative_" : "")$(get(m1, "case", "runs"))_N$(get(m1, "N", "?"))_compare.png")
    savefig(plot(p; size = (520, 420), margin = 4Plots.mm, left_margin = 8Plots.mm), out)
    println("\nwrote ", out)
end

main(ARGS)
