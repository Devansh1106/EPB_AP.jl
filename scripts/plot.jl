# Plot CSV files written by EPB_AP (solution, snapshot or history files). Does not run the solver.
#
#   julia --project=scripts scripts/plot.jl data/<file>.csv [more.csv ...] [out.png] [--var rho,u,phi] [--zoom x1,x2[,y1,y2]]
#         [--inset topright|top|bottom|center|...] [--legend bottomleft|...] [--no-initial] [--no-exact]
#
# The first column is the abscissa; every other column gets a panel (--var: only these columns, in this order), its
# axis labelled by the variable's symbol (ρ, u, ϕ, ...). The references are black lines of width 1: the exact solution (column `<name>_exact`) solid, the
# initial data (`<name>_init`) solid as well when the panel has no exact solution, dashed otherwise. The numerical
# solutions each get their own colour and marker (scripts/style.jl), labelled by scheme (FOPD, FODD, SODD); --no-initial
# (or --no-init), --no-exact omit those columns. --legend places the legend (default best, bottomleft with --zoom).
# --zoom x1,x2 adds to each panel an inset magnifying [x1, x2] (y-range fitted; x1,x2,y1,y2 sets it, the same in all
# panels), the region marked by a dotted box; --inset places it (top right by default).
# Several files are drawn in the same panels, each curve labelled by its scheme and the metadata that differs: case, N,
# λ, t; if that is not enough, any other header entry that differs (limiter, θ, timestep, ...); else the file name. There is
# no title. Initial data are drawn once per case, exact solutions once per case and time. Default output: the first file with .png (`_compare.png` for several).
# The exact solution of riemann (and riemann_speed) is the λ → 0 one: plotting it for λ > 1e-4 prints a warning.
using Plots
include("style.jl")

function read_data(path)
    meta, names, rows = Pair{String,String}[], String[], Vector{Float64}[]
    for line in eachline(path)
        if startswith(line, "#")
            k, v = strip.(split(line[2:end], "="; limit = 2))
            push!(meta, k => v)
        elseif isempty(names)
            names = String.(split(line, ","))
        elseif !isempty(line)
            push!(rows, parse.(Float64, split(line, ",")))
        end
    end
    return Dict(meta), names, reduce(vcat, permutedims.(rows))
end

const KEYS = ("case", "scheme", "N", "lambda", "t")
const SHOW = Dict("lambda" => "λ")
entry(meta, k) = k == "scheme" ? scheme_label(get(meta, k, "?")) : "$(get(SHOW, k, k))=$(get(meta, k, "?"))"
tag(meta, keys) = join([entry(meta, k) for k in keys], ", ")

# Legend label of each file: its scheme (FOPD, FODD, SODD) and the other KEYS that differ between the files; while two
# labels are the same, the other header entries that differ too (but the run's step count and status); if that is
# not enough, the file name.
function file_labels(paths, files)
    differs(k) = length(unique(get(m, k, "?") for (m, _, _) in files)) > 1
    keys = filter(k -> differs(k) && k != "scheme", collect(KEYS))
    label(m, keys) = join(filter(!isempty, [haskey(m, "scheme") ? entry(m, "scheme") : ""; tag(m, keys)]), ", ")
    labs = [label(m, keys) for (m, _, _) in files]
    if !allunique(labs)
        more = sort(filter(k -> differs(k) && !(k in KEYS) && !(k in ("steps", "status")),
                           unique(k for (m, _, _) in files for k in Base.keys(m))))
        labs = [label(m, [keys; more]) for (m, _, _) in files]
    end
    return allunique(labs) && all(!isempty, labs) ? labs : basename.(paths)
end

# Cases whose "exact" solution is the quasineutral (λ → 0) one, and the largest λ at which it is still a fair reference.
const LIMIT_EXACT = ("riemann", "riemann_speed")
const LIMIT_LAMBDA = 1e-4

function warn_limit_exact(paths, files)
    for (path, (m, names, _)) in zip(paths, files)
        get(m, "case", "") in LIMIT_EXACT && any(endswith("_exact"), names) || continue
        λ = tryparse(Float64, get(m, "lambda", ""))
        λ !== nothing && λ > LIMIT_LAMBDA &&
            @warn "$(basename(path)): the exact solution of $(m["case"]) is the quasineutral (λ → 0) one, " *
                  "[1, (7.8)–(7.10)]; at λ = $λ > $LIMIT_LAMBDA it is not the solution of this problem " *
                  "(--no-exact omits it)"
    end
end

# Draws the panel of column `c`: the numerical solutions, then on top the references (initial data, dashed if an exact
# solution follows; exact). With `window = (x1, x2)`, only the points in [x1, x2] (and one on each side). `kw` go to
# every series (e.g. the subplot of an inset). Returns the values drawn, for the axis limits.
function draw!(p, c, files, labels; init, exact, window = nothing, kw...)
    drawn, ys, k = Set{Any}(), Float64[], 0
    function sel(x)
        window === nothing && return eachindex(x)
        i = findall(v -> window[1] <= v <= window[2], x)
        return isempty(i) ? i : max(first(i) - 1, 1):min(last(i) + 1, length(x))
    end
    extra(suffix, m) = (suffix, suffix == "_init" ? (m["case"],) : (m["case"], get(m, "t", "")))
    has_exact = exact && any(((_, names, _),) -> c * "_exact" in names, files)
    function others!(suffix, name)
        for (i, (m, names, data)) in enumerate(files)
            j = findfirst(==(c * suffix), names)
            (j === nothing || extra(suffix, m) in drawn) && continue
            push!(drawn, extra(suffix, m))
            lab = count(d -> d[1] == suffix, drawn) == 1 ? name : "$name, " * labels[i]
            r = sel(data[:, 1])
            reference!(p, data[r, 1], data[r, j]; label = lab, ls = suffix == "_init" && has_exact ? :dash : :solid, kw...)
            (suffix == "_init" || window !== nothing) && append!(ys, data[r, j])
        end
    end
    for (i, (m, names, data)) in enumerate(files)
        j = findfirst(==(c), names)
        j === nothing && continue
        r = sel(data[:, 1])
        append!(ys, data[r, j])
        curve!(p, data[r, 1], data[r, j], k += 1; label = labels[i], kw...)
    end
    init && others!("_init", "Initial")
    exact && others!("_exact", "Exact")
    return ys
end

# Inset positions (in the panel's plot area) for --inset.
const INSETS = Dict("bottom" => (0.0, 0.06, :bottom, :hcenter), "bottomleft" => (0.06, 0.06, :bottom, :left),
                    "bottomright" => (0.06, 0.06, :bottom, :right), "top" => (0.0, 0.06, :top, :hcenter),
                    "topleft" => (0.06, 0.06, :top, :left), "topright" => (0.06, 0.06, :top, :right),
                    "center" => (0.0, 0.0, :vcenter, :hcenter))

function plot_files(paths; out = replace(paths[1], r"\.csv$" => length(paths) == 1 ? ".png" : "_compare.png"),
                    vars = nothing, init = true, exact = true, zoom = nothing, inset = "topright",
                    legend = zoom === nothing ? :best : :bottomleft)
    files = read_data.(paths)
    labels = file_labels(paths, files)
    exact && warn_limit_exact(paths, files)
    _, names1, _ = files[1]
    cols = filter(n -> n != names1[1] && !endswith(n, "_exact") && !endswith(n, "_init"), names1)
    if vars !== nothing
        bad = setdiff(vars, cols)
        isempty(bad) || error("unknown variable(s) $(join(bad, ", ")) in --var; available: $(join(cols, ", "))")
        cols = vars
    end
    panels = map(cols) do c
        p = plot(; xlabel = names1[1], ylabel = symbol(c), legend)
        ys = draw!(p, c, files, labels; init, exact)
        c in ("dt", "ap_residual") && all(>(0), filter(isfinite, ys)) && plot!(p; yscale = :log10)
        lo, hi = extrema(filter(isfinite, ys))
        w = max(1, abs(hi))
        hi - lo < 1e-9w && ylims!(p, (lo - 1e-6w, hi + 1e-6w))   # constant up to round-off (e.g. mass)
        p
    end
    nc = min(length(panels), 3); nr = cld(length(panels), nc)   # near-square panels, at most 3 per row
    fig = plot(panels...; layout = (nr, nc), size = (420nc, 400nr), margin = 4Plots.mm, left_margin = 10Plots.mm, bottom_margin = 8Plots.mm)
    if zoom !== nothing
        # One inset per panel: [x1, x2] magnified, its y-range [y1, y2] or fitted to the curves there; a dotted box
        # marks the region.
        x0, y0, a1, a2 = INSETS[inset]
        for (i, c) in enumerate(cols)
            ys = filter(isfinite, draw!(plot(), c, files, labels; init, exact, window = zoom))
            isempty(ys) && error("no data of $c in --zoom $(zoom[1]),$(zoom[2])")
            lo, hi = length(zoom) == 4 ? zoom[3:4] : extrema(ys)
            pad = length(zoom) == 4 ? 0.0 : max(0.05(hi - lo), 1e-12max(1, abs(hi)))
            plot!(fig, [zoom[1], zoom[2], zoom[2], zoom[1], zoom[1]], [lo, lo, hi, hi, lo] .+ pad .* [-1, -1, 1, 1, -1];
                  subplot = i, color = :gray, ls = :dot, lw = 1, label = "")
            j = length(cols) + i
            plot!(fig, Float64[], Float64[]; inset = (i, bbox(x0, y0, 0.4, 0.32, a1, a2)), subplot = j, label = "")
            draw!(fig, c, files, labels; init, exact, window = zoom, subplot = j)
            plot!(fig; subplot = j, legend = false, xlims = zoom[1:2], ylims = (lo - pad, hi + pad), tickfontsize = 6,
                  framestyle = :box, background_color_inside = :white)
        end
    end
    savefig(fig, out)
    println("wrote ", out)
end

const FLAGS = ("--no-init", "--no-initial", "--no-exact")
const LEGENDS = ("best", "topleft", "top", "topright", "left", "right", "bottomleft", "bottom", "bottomright",
                 "outertopright", "outerright", "outerbottom", "none")
isout(a) = occursin(r"\.(png|pdf|svg)$", a)
usage = "usage: julia --project=scripts scripts/plot.jl data/<file>.csv [more.csv ...] [out.png] [--var rho,u,phi] " *
        "[--zoom x1,x2[,y1,y2]] [--inset topright|top|bottom|center|...] [--legend bottomleft|...] [--no-init(ial)] [--no-exact]"
argv, vals = copy(ARGS), Dict{String,String}()
for k in ("--var", "--zoom", "--inset", "--legend")
    i = findfirst(==(k), argv)
    i === nothing && continue
    (i == length(argv) || startswith(argv[i+1], "--")) && error("$k needs a value\n" * usage)
    vals[k] = argv[i+1]
    deleteat!(argv, i:i+1)
end
vars = haskey(vals, "--var") ? String.(filter(!isempty, strip.(split(vals["--var"], ",")))) : nothing
zoom = haskey(vals, "--zoom") ? Tuple(parse.(Float64, split(vals["--zoom"], ","))) : nothing
zoom === nothing || (length(zoom) in (2, 4) && zoom[1] < zoom[2] && (length(zoom) == 2 || zoom[3] < zoom[4])) ||
    error("--zoom takes x1,x2 (x1 < x2) or x1,x2,y1,y2 (y1 < y2: one y-range for all panels)\n" * usage)
inset = get(vals, "--inset", "topright")
haskey(INSETS, inset) || error("--inset is one of $(join(sort(collect(keys(INSETS))), ", "))\n" * usage)
legend = get(vals, "--legend", nothing)
legend === nothing || legend in LEGENDS || error("--legend is one of $(join(LEGENDS, ", "))\n" * usage)
flags, args = filter(startswith("--"), argv), filter(!startswith("--"), argv)
inputs, outs = filter(!isout, args), filter(isout, args)
(isempty(inputs) || length(outs) > 1) && error(usage)
for f in flags
    f in FLAGS || error("unknown option $f\n" * usage)
end
plot_files(inputs; (isempty(outs) ? () : (out = outs[1],))..., (legend === nothing ? () : (legend = Symbol(legend),))...,
           vars, zoom, inset, init = isempty(intersect(flags, ("--no-init", "--no-initial"))), exact = !("--no-exact" in flags))
