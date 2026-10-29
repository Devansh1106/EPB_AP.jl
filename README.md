# EPB_AP.jl

[![Build Status](https://github.com/Devansh1106/EPB_AP.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Devansh1106/EPB_AP.jl/actions/workflows/CI.yml?query=branch%3Amain)

Collocated finite volume schemes for the Euler–Poisson–Boltzmann system (Arun & Ghorai): the first-order
shift scheme of §3.3 (`Scheme()`), the density-dissipation scheme of §3.4 (`Scheme(kind = :dissipation)`) and
the second-order MUSCL/IMEX scheme of §4 (`SecondOrder()`). Equations, test cases and choices are listed in
[SPEC.md](SPEC.md).

## Install

```julia
pkg> add https://github.com/Devansh1106/EPB_AP.jl
```
or, from a clone: `julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'`.

## Run a test case

```julia
using EPB_AP
run_case(:step)                                      # writes data/step_dissipation_N200_lambda0.1_T1.0*.csv
run_case(:smooth; lambda = 1e-2, scheme = Scheme())  # override any default
run_case(:smooth; scheme = SecondOrder())            # second order (limiter = :minmod or :none, c_D, cfl)
```
In a parameter file, second order is `kind = "second"` under `[scheme]`; keys of the other order (`theta`, `kappa_*`
or `limiter`, `c_D`, `cfl`) are ignored, so only `kind` needs to change.
`run_case` prints the files it writes and warns if the run stopped before `T` (a collapse); the files then hold the
state at the time reached (`# t` in the header), while their names keep `T`.
Cases: `smooth, steep, nearvacuum, bump, step, expansion, soliton, ap, simplewave, riemann, riemann_speed, fivebranch, sevenbranch, shocktube, revshocktube, positivity, plasmaexpansion`.

## Parameters

Each case reads its parameters (domain, boundary, N, λ, T, scheme, time step) from `params/<case>.toml`.
Edit that file, or copy it and run the copy: `run_case("my_step.toml")`. Keywords override the file.

## Plot

All scripts, with every option: [scripts/README.md](scripts/README.md).

```sh
julia --project=scripts -e 'using Pkg; Pkg.instantiate()'    # once
julia --project=scripts scripts/plot.jl data/step_dissipation_N200_lambda0.1_T1.0.csv
```

Each panel shows the numerical solution, the initial data and, if known, the exact solution. Plot style (all scripts,
`scripts/style.jl`): references are black lines of width 1 (exact solution and reference slopes solid; initial data solid, or dashed
when an exact solution is drawn too); every other curve
has its own colour and marker, width 1.2; axes use the symbols ρ, u, ϕ; no titles.
Several files are drawn in the same panels, e.g. two schemes or grids, with an optional output name:
`scripts/plot.jl data/a.csv data/b.csv compare.png`. `--no-initial` (or `--no-init`) and `--no-exact` leave out the initial data and the exact solution; `--var rho,phi` plots only those
columns, in that order (default: all); `--zoom x1,x2` adds to each panel an inset magnifying [x1, x2] (`x1,x2,y1,y2`
also fixes its y-range), placed by `--inset topright|top|bottom|center|...` (default top right); `--legend bottomleft|...` moves the
legend (default: best place, bottom left with `--zoom`). Curves are labelled by scheme: FOPD, FODD, SODD. Legends name the schemes FOPD (shift), FODD (dissipation) and SODD (second).
The exact solution of `riemann` (and `riemann_speed`) is the quasineutral (λ → 0) one, so plotting it for λ > 1e-4 prints a warning.

### Plot with matplotlib (for a manuscript)

`scripts/plot.py` draws the same plots with matplotlib (needs `numpy`, `matplotlib`), with the same options, sized
for LaTeX: Computer Modern rendered by LaTeX (`--no-usetex`: matplotlib's own; LaTeX needs the `type1cm` package),
9 pt labels, 8 pt ticks and legend, full frame with ticks inside, no grid. `--width` is the figure width in inches
(default 3.4 for one panel, 6.5 for more), so include it at that size; output PNG at 300 dpi (or name `.pdf`, `.svg`).
With `--zoom`, a grey box marks the region and lines join it to the inset's corners.
```sh
python3 scripts/plot.py data/riemann_shift_N200_lambda0.0001_T50.0.csv data/riemann_second_N200_lambda0.0001_T50.0.csv \
        riemann_rho.png --var rho --zoom 45,62,0.69,0.73 --no-initial
```

## Convergence

```sh
julia --project=scripts scripts/convergence.jl simplewave --scheme shift --N 100,200,400,800
julia --project=scripts scripts/convergence.jl soliton --T 2 --dtmax 0.4 --var rho,m   # only some variables
julia --project=scripts scripts/convergence.jl smooth --ref 3200     # no exact solution: fine-grid reference
julia --project=scripts scripts/convergence.jl simplewave --scheme second --fixed 0.4   # Δt = 0.4Δx
julia --project=scripts scripts/convergence.jl smooth --scheme second --limiter none --c_D 0 --ref 3200
julia --project=scripts scripts/convergence.jl smooth --scheme second --plot yes --slope 1,2   # log–log plot
```
`--scheme` keeps the file's settings when the kind is unchanged; unknown options are an error.
`--plot yes` writes `convergence_<case>_<scheme>[_<limiter>]_lambda<λ>.png` (or `--plot out.png`): L1, L2 and L∞
errors against N per variable, with black reference slopes (default: the scheme's order; with several, the second
and later ones dashed, dotted).

## Energy stability

```sh
julia --project=scripts scripts/energy.jl smooth --plot yes      # shift, dissipation, second-none, second-minmod; λ = 1 … 1e-6
julia --project=scripts scripts/energy.jl ap --schemes dissipation,second --lambda 1,1e-4 --N 400
```
For each scheme and λ: E(T) − E(0), the largest step change max(Eⁿ⁺¹ − Eⁿ) of the energy (14) and the number of steps
where it increases. `--plot yes` writes `energy_<case>_N<N>.png`: E(t) − E(0) against t, one panel per scheme.
`--relative` uses E(t)/E(0) instead (the energy (14) unshifted, so with E(0) < 0 a decreasing energy gives a ratio
above 1). Given history files of runs instead of a case, nothing is run: one plot of those runs and, if the case has
one, the exact solution's energy (black), e.g. the soliton over two periods:
```sh
julia --project=. -e 'using EPB_AP; T = 100/1.3; run_case(:soliton; N = 2000, T); run_case(:soliton; N = 2000, T, scheme = SecondOrder())'
julia --project=scripts scripts/energy.jl data/soliton_dissipation_N2000_lambda1.0_T76.92307692307692_history.csv \
      data/soliton_second_N2000_lambda1.0_T76.92307692307692_history.csv --relative --plot soliton_energy.png
```
