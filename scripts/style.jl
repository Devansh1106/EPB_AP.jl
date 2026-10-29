# Plot style shared by the scripts: a reference (exact solution, reference slope) is a solid black line of width 1;
# every other curve has its own colour and marker and width 1.2. No titles.

const MARKERS =(:circle, :rect, :utriangle, :diamond, :dtriangle, :star5, :pentagon, :hexagon, :xcross, :cross)
const REF_LW = 1.0
const LW = 1.2

"Axis symbol of a solution variable (CSV column or convergence variable)."
const SYMBOLS = Dict("rho" => "ρ", "u" => "u", "phi" => "ϕ", "m" => "ρu", "energy" => "E", "min_rho" => "min ρ",
                     "mass" => "mass", "ap_residual" => "‖e^ϕ − ρ‖∞", "dt" => "Δt", "cfl_ratio" => "CFL ratio")
symbol(name) = get(SYMBOLS, name, name)

# Legend name of a scheme (`scheme` header entry): shift (§3.3), dissipation (§3.4), second order (§4).
const SCHEME_LABELS = Dict("shift" => "FOPD", "dissipation" => "FODD", "second" => "SODD")
scheme_label(name) = get(SCHEME_LABELS, name, name)

reference!(p, x, y; label = "reference", kw...) = plot!(p, x, y; color = :black, lw = REF_LW, ls = :solid, label, kw...)

"""
    curve!(p, x, y, k; label, nmarks = 15)

Curve number `k` (its colour and marker): the line through all points, markers on at most `nmarks` of them so that
fine grids stay readable (staggered with `k`, so that curves on the same grid do not hide each other's markers),
and one legend entry showing both.
"""
function curve!(p, x, y, k; label, nmarks = 15, kw...)
    c, m = k, MARKERS[mod1(k, length(MARKERS))]
    plot!(p, x, y; color = c, lw = LW, label = "", kw...)
    n = length(x)
    step = (n - 1) / (nmarks - 1)
    idx = n <= nmarks ? eachindex(x) : unique(clamp.(round.(Int, 1 .+ mod((k - 1) * 0.3step, step) .+ step .* (0:nmarks-1)), 1, n))
    idx = filter(i -> isfinite(y[i]), idx)
    scatter!(p, x[idx], y[idx]; color = c, marker = m, markersize = 4, markerstrokewidth = 0.5, label = "", kw...)
    plot!(p, [NaN], [NaN]; color = c, lw = LW, marker = m, markersize = 4, markerstrokewidth = 0.5, label, kw...)
end
