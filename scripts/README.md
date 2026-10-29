# scripts/

Scripts for plotting runs and for convergence and energy studies. Run them from the repository root.

| Script | Language | Runs the solver | Purpose |
|---|---|---|---|
| [`plot.jl`](#plotjl--plotpy) | Julia (Plots.jl, GR) | no | plot CSV files written by `run_case` |
| [`plot.py`](#plotjl--plotpy) | Python (matplotlib) | no | the same plots, styled for a LaTeX manuscript |
| [`convergence.jl`](#convergencejl) | Julia | yes | errors and orders on refined grids |
| [`energy.jl`](#energyjl) | Julia | yes (or reads history files) | the discrete energy (14) along runs |
| `style.jl`, `style.py` | | | plot style shared by the Julia / Python scripts |

## Setup

```sh
julia --project=scripts -e 'using Pkg; Pkg.instantiate()'   # once: EPB_AP (from ..) and Plots
pip install numpy matplotlib                                 # for plot.py
```
`plot.py` renders text with LaTeX when `latex` is installed (it needs the `type1cm` package, e.g.
`sudo dnf install texlive-type1cm`); otherwise, or with `--no-usetex`, matplotlib's own Computer Modern.

## Input files

`run_case` writes to `data/`:
- `<case>_<scheme>_N<N>_lambda<λ>_T<T>.csv`: the solution at the final time. Header lines `# key = value` (case,
  scheme, N, domain, bc, lambda, T, t reached, steps, status, timestep, scheme settings), then the columns
  `x, rho, u, phi`, the initial data `rho_init, u_init, phi_init` and, if the case has one, the exact solution
  `rho_exact, u_exact, phi_exact`.
- `..._history.csv`: one row per time step, columns `t, dt, energy, min_rho, mass, ap_residual, cfl_ratio`.
- snapshot files, if the case asks for them, in the solution format.

## Plot style

Both `style.jl` and `style.py`:
- References (the exact solution, the initial data, reference slopes) are black lines of width 1. The initial data
  are dashed when an exact solution is drawn as well.
- Every other curve has its own colour and marker and width 1.2, with markers on about 15 points, staggered between
  curves so that they do not hide each other.
- Axes are labelled by the symbols ρ, u, ϕ, ρu, E, min ρ, Δt, ‖e^ϕ − ρ‖∞. There are no titles.
- Legends name the schemes FOPD (shift, §3.3), FODD (dissipation, §3.4) and SODD (second order, §4).

`style.py` also sets the manuscript look: Computer Modern, 9 pt labels, 8 pt ticks and legend, a full frame with
ticks inside on all sides, no grid, a thin black legend frame and 300 dpi PNG.

## plot.jl / plot.py

Plot one or more CSV files (solution, snapshot or history files) in the same panels.

```sh
julia --project=scripts scripts/plot.jl data/<file>.csv [more.csv ...] [out.png] [options]
python3 scripts/plot.py data/<file>.csv [more.csv ...] [out.png] [options]
```

The first column is the abscissa, and every other column gets a panel, at most 3 per row. Each curve is labelled by
its scheme plus whatever differs between the files: first case, N, λ, t; then any other header entry (limiter, θ,
timestep, ...); if that is still not enough, the file name. Initial data are drawn once per case, exact solutions once
per case and time.

| Option | Meaning |
|---|---|
| `out.png` | output file. Default: the first file with `.png`, or `_compare.png` for several files. `plot.jl` also accepts `.pdf`, `.svg`; `plot.py` also `.pgf`, `.eps` |
| `--var rho,u,phi` | plot only these columns, in this order (default: all) |
| `--zoom x1,x2` | add to each panel an inset magnifying [x1, x2], its y-range fitted to the curves there |
| `--zoom x1,x2,y1,y2` | the same, with the inset's y-range fixed (the same in all panels) |
| `--inset POS` | inset position: `topright` (default), `top`, `topleft`, `bottom`, `bottomleft`, `bottomright`, `center` |
| `--legend POS` | legend position: `best` (default), `bottomleft` (default with `--zoom`), `topleft`, `top`, `topright`, `left`, `right`, `bottom`, `bottomright`, `outertopright`, `outerright`, `outerbottom`, `none` |
| `--no-initial` (`--no-init`) | leave out the initial data |
| `--no-exact` | leave out the exact solution |
| `--width W` | (`plot.py` only) figure width in inches. Default: 3.4 (one column) for one panel, 6.5 (`\textwidth`) for more. Include the figure at this width so that its fonts match the text |
| `--no-usetex` | (`plot.py` only) use matplotlib's mathtext instead of LaTeX |

Notes:
- The zoomed region is marked in the panel: `plot.jl` draws a dotted box; `plot.py` draws a grey box with lines to
  the inset's corners.
- The exact solution of `riemann` and `riemann_speed` is the quasineutral (λ → 0) one, [1, (7.8)–(7.10)]. Plotting
  it for λ > 1e-4 prints a warning.
- Δt and the AP residual are drawn on a log axis: always in `plot.jl`; in `plot.py` only when they span a decade
  or more. Nearly constant columns (e.g. mass) get a y-range around their value.
- `plot.py` takes about 4 s; `plot.jl` about 35 s (loading Plots).

Examples:
```sh
# FOPD vs SODD for riemann, ρ only, inset on the overshoot behind the shock
python3 scripts/plot.py data/riemann_shift_N200_lambda0.0001_T50.0.csv data/riemann_second_N200_lambda0.0001_T50.0.csv \
        riemann_rho.png --var rho --zoom 45,62,0.69,0.73 --no-initial
# the same for u and ϕ
...  riemann_u.png --var u --zoom 45,62,0.33,0.37 --inset bottom --legend topleft --no-initial
...  riemann_phi.png --var phi --zoom 45,62,-0.36,-0.32 --no-initial
# a history file: one panel per column
julia --project=scripts scripts/plot.jl data/soliton_second_N2000_lambda1.0_T76.92307692307692_history.csv
```

## convergence.jl

Errors of one test case on a sequence of refined grids, and the observed orders.

```sh
julia --project=scripts scripts/convergence.jl <case | file.toml> [options]
```

| Option | Meaning |
|---|---|
| `--scheme shift\|dissipation\|second` | scheme kind (default: the case file's; the file's settings are kept when the kind is the same) |
| `--limiter minmod\|none` | limiter of the second-order scheme |
| `--c_D c` | dissipation coefficient of the second-order scheme |
| `--lambda λ` | λ (default: the file's) |
| `--T T` | final time (default: the file's) |
| `--N 100,200,400,800` | grids (this is the default) |
| `--dtmax c` | cap the scheme's Δt rule at c·Δx |
| `--fixed c` | fixed Δt = c·Δx (the §5.4 order studies use 0.4) |
| `--var rho,u,phi,m` | variables, m = ρu (default: all four) |
| `--ref Nref` | grid of the reference solution when the case has no exact solution. Default 4·max N; it must be a multiple of every N |
| `--plot yes\|out.png` | log–log plot of the L1, L2, L∞ errors against N, one panel per variable. `yes` writes `convergence_<case>_<scheme>[_<limiter>]_lambda<λ>.png` |
| `--slope 1,2` | reference slopes in the plot (default: the scheme's order). The first is solid, the next dashed, dotted. Needs `--plot` |

Errors are measured against the exact solution when the case has one. Otherwise the reference is the same scheme on
Nref cells, restricted to each grid by cell averaging. The norms are normalized: (Δx/|Ω|)Σ|e|, ((Δx/|Ω|)Σe²)^½ and
max|e|, so L1 ≤ L2 ≤ L∞. One table is printed per variable, and each grid is solved once. Unknown options are an
error.

Examples:
```sh
julia --project=scripts scripts/convergence.jl simplewave --scheme shift --N 100,200,400,800
julia --project=scripts scripts/convergence.jl soliton --T 2 --dtmax 0.4 --var rho,m
julia --project=scripts scripts/convergence.jl smooth --scheme second --limiter none --c_D 0 --ref 3200
julia --project=scripts scripts/convergence.jl smooth --scheme second --plot yes --slope 1,2
```

## energy.jl

The discrete energy (14), E = Σ ρu²/2 + e^ϕ(ϕ − 1) + λ²/2 |∂ϕ|², along runs. It works in two modes.

**Case mode** runs a case for several schemes and λ:
```sh
julia --project=scripts scripts/energy.jl <case | file.toml> [options]
```

| Option | Meaning |
|---|---|
| `--schemes shift,dissipation,second-none,second-minmod` | schemes (this is the default); `kind[-limiter]`, keeping the file's settings when the kind is the file's |
| `--lambda 1,1e-2,1e-4,1e-6` | values of λ (this is the default) |
| `--N N` | grid (default: the file's) |
| `--T T` | final time (default: the file's) |
| `--plot yes\|out.png` | E(t) − E(0) against t, one panel per scheme, one curve per λ. `yes` writes `energy_<case>_N<N>.png` |
| `--relative` | E(t)/E(0) instead of E(t) − E(0), in the table and the plot |

**File mode** reads `_history.csv` files written by `run_case` and runs nothing:
```sh
julia --project=scripts scripts/energy.jl <run>_history.csv [more_history.csv ...] [--plot out.png] [--relative]
```
It draws one plot of all the runs (legend FOPD, FODD, SODD). If the case has an exact solution, its energy on the
first file's grid is the black reference "Exact"; it is left out, with a warning, when its E(0) differs from the
run's by more than 1%. Default output: `energy_[relative_]<case>_N<N>_compare.png`.

For each run the table gives:
- E(T) − E(0), or E(T)/E(0) with `--relative`;
- the largest step change max(Eⁿ⁺¹ − Eⁿ);
- the number of steps where the energy increases beyond round-off (1e3·eps·max(1, |E⁰|)).

Energy stability means no increasing step: Theorems 3.5 and 3.16 for first order. For second order no fully discrete
inequality is claimed (§4.13), so the table is an observation.

`--relative` divides by E(0), which is negative for e.g. the soliton (E(0) ≈ −47.5). A decreasing energy then gives a
ratio above 1.

Examples:
```sh
julia --project=scripts scripts/energy.jl smooth --plot yes
julia --project=scripts scripts/energy.jl ap --schemes dissipation,second --lambda 1,1e-4 --N 400
# soliton over two periods, FODD and SODD
julia --project=. -e 'using EPB_AP; T = 100/1.3; run_case(:soliton; N = 2000, T); run_case(:soliton; N = 2000, T, scheme = SecondOrder())'
julia --project=scripts scripts/energy.jl data/soliton_dissipation_N2000_lambda1.0_T76.92307692307692_history.csv \
      data/soliton_second_N2000_lambda1.0_T76.92307692307692_history.csv --relative --plot soliton_energy.png
```
