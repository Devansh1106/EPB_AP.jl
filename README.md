# EPB_AP.jl

[![Build Status](https://github.com/Devansh1106/EPB_AP.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Devansh1106/EPB_AP.jl/actions/workflows/CI.yml?query=branch%3Amain)

First-order collocated finite volume schemes for the Euler–Poisson–Boltzmann system
(Arun & Ghorai): the shift scheme of §3.3 (`Scheme()`) and the density-dissipation scheme of
§3.4 (`Scheme(kind = :dissipation)`). Equations, test cases and choices are listed in [SPEC.md](SPEC.md).

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
```
Cases: `smooth, steep, nearvacuum, bump, step, expansion, soliton, ap, simplewave, riemann, riemann_speed`.

## Parameters

Each case reads its parameters (domain, boundary, N, λ, T, scheme, time step) from `params/<case>.toml`.
Edit that file, or copy it and run the copy: `run_case("my_step.toml")`. Keywords override the file.

## Plot

```sh
julia --project=scripts -e 'using Pkg; Pkg.instantiate()'    # once
julia --project=scripts scripts/plot.jl data/step_dissipation_N200_lambda0.1_T1.0.csv
```

Each panel shows the initial data (dotted), the numerical solution and, if known, the exact solution (dashed).

## Convergence

```sh
julia --project=scripts scripts/convergence.jl simplewave --scheme shift --N 100,200,400,800
julia --project=scripts scripts/convergence.jl soliton --T 2 --dtmax 0.4 --var rho,m   # only some variables
julia --project=scripts scripts/convergence.jl smooth --dtmax 0.4 --ref 3200   # no exact solution: fine-grid reference
```
