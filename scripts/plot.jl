# Plot any CSV file written by EPB_AP (solution, snapshot or history file). Does not run the solver.
#
#   julia --project=scripts scripts/plot.jl data/<file>.csv [out.png]
#
# The first column is the abscissa; every other column gets a panel, and a column `<name>_exact`
# is drawn dashed in the panel of `<name>`. The header metadata becomes the title.
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

function plot_file(path; out = replace(path, r"\.csv$" => ".png"))
    meta, names, data = read_data(path)
    x = data[:, 1]
    cols = filter(n -> n != names[1] && !endswith(n, "_exact"), names)
    panels = map(cols) do c
        y = data[:, findfirst(==(c), names)]
        p = plot(x, y; label = "numerical", xlabel = names[1], ylabel = c, lw = 1.5,
                 yscale = c in ("dt", "ap_residual") && all(>(0), filter(isfinite, y)) ? :log10 : :identity)
        lo, hi = extrema(filter(isfinite, y))
        w = max(1, abs(hi))
        hi - lo < 1e-9w && ylims!(p, (lo - 1e-6w, hi + 1e-6w))   # constant up to round-off (e.g. mass)
        e = findfirst(==(c * "_exact"), names)
        e === nothing || plot!(p, x, data[:, e]; label = "exact", ls = :dash, lw = 1.5)
        p
    end
    title = join(filter(!isempty, [get(meta, k, "") for k in ("case", "scheme")]), ", ") *
            "  N=$(get(meta, "N", "?")), ε=$(get(meta, "eps", "?")), t=$(get(meta, "t", "?"))"
    fig = plot(panels...; layout = length(panels), size = (450 * min(length(panels), 3),
               330 * cld(length(panels), 3)), plot_title = title, plot_titlefontsize = 10, margin = 4Plots.mm, left_margin = 10Plots.mm, bottom_margin = 8Plots.mm)
    savefig(fig, out)
    println("wrote ", out)
end

isempty(ARGS) && error("usage: julia --project=scripts scripts/plot.jl data/<file>.csv [out.png]")
plot_file(ARGS[1]; (length(ARGS) > 1 ? (out = ARGS[2],) : ())...)
