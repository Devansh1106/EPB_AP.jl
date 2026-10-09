# Plot CSV files written by EPB_AP (solution, snapshot or history files). Does not run the solver.
#
#   julia --project=scripts scripts/plot.jl data/<file>.csv [more.csv ...] [out.png] [--no-init] [--no-exact]
#
# The first column is the abscissa; every other column gets a panel. In the panel of `<name>`, a column
# `<name>_init` (initial data) is drawn dotted and a column `<name>_exact` dashed; --no-init, --no-exact omit them.
# Several files are drawn in the same panels: metadata shared by all files becomes the title, metadata
# that differs (scheme, N, λ, t, ...) labels each curve. Initial data are drawn once per case, exact
# solutions once per case and time. Default output: the first file with .png (`_compare.png` for several).
using Plots

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
tag(meta, keys) = join(["$(get(SHOW, k, k))=$(get(meta, k, "?"))" for k in keys], ", ")

function plot_files(paths; out = replace(paths[1], r"\.csv$" => length(paths) == 1 ? ".png" : "_compare.png"),
                    init = true, exact = true)
    files = read_data.(paths)
    common = [k for k in KEYS if length(unique(get(m, k, "?") for (m, _, _) in files)) == 1]
    differ = [k for k in KEYS if !(k in common)]
    _, names1, _ = files[1]
    cols = filter(n -> n != names1[1] && !endswith(n, "_exact") && !endswith(n, "_init"), names1)
    panels = map(cols) do c
        p = plot(; xlabel = names1[1], ylabel = c)
        drawn, ys = Set{Any}(), Float64[]
        for (m, names, data) in files                      # initial data and exact solutions first
            x = data[:, 1]
            for (suffix, key, style) in (("_init", (m["case"],), (label = "initial", ls = :dot, color = :gray)),
                                         ("_exact", (m["case"], get(m, "t", "")), (label = "exact", ls = :dash, color = :red)))
                j = findfirst(==(c * suffix), names)
                (j === nothing || (suffix, key) in drawn || !(suffix == "_init" ? init : exact)) && continue
                push!(drawn, (suffix, key))
                lab = count(d -> d[1] == suffix, drawn) == 1 ? style.label : "$(style.label), " * tag(m, differ)
                plot!(p, x, data[:, j]; lw = 1.5, style..., label = lab)
                suffix == "_init" && append!(ys, data[:, j])
            end
        end
        for (k, (m, names, data)) in enumerate(files)
            j = findfirst(==(c), names)
            j === nothing && continue
            y = data[:, j]
            append!(ys, y)
            lab = length(files) == 1 ? "numerical" : tag(m, differ)
            plot!(p, data[:, 1], y; label = lab, lw = 1.5, color = length(files) == 1 ? :blue : k,
                  yscale = c in ("dt", "ap_residual") && all(>(0), filter(isfinite, y)) ? :log10 : :identity)
        end
        lo, hi = extrema(filter(isfinite, ys))
        w = max(1, abs(hi))
        hi - lo < 1e-9w && ylims!(p, (lo - 1e-6w, hi + 1e-6w))   # constant up to round-off (e.g. mass)
        p
    end
    title = tag(files[1][1], common)
    nc = min(length(panels), 3); nr = cld(length(panels), nc)   # near-square panels, at most 3 per row
    fig = plot(panels...; layout = (nr, nc), size = (420nc, 400nr + 30), plot_title = title, plot_titlefontsize = 10, margin = 4Plots.mm, left_margin = 10Plots.mm, bottom_margin = 8Plots.mm)
    savefig(fig, out)
    println("wrote ", out)
end

const FLAGS = ("--no-init", "--no-exact")
isout(a) = occursin(r"\.(png|pdf|svg)$", a)
flags, args = filter(startswith("--"), ARGS), filter(!startswith("--"), ARGS)
inputs, outs = filter(!isout, args), filter(isout, args)
usage = "usage: julia --project=scripts scripts/plot.jl data/<file>.csv [more.csv ...] [out.png] [--no-init] [--no-exact]"
(isempty(inputs) || length(outs) > 1) && error(usage)
for f in flags
    f in FLAGS || error("unknown option $f\n" * usage)
end
plot_files(inputs; (isempty(outs) ? () : (out = outs[1],))..., init = !("--no-init" in flags), exact = !("--no-exact" in flags))
