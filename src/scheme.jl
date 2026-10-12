"""
    Scheme(; kind = :shift, theta = 0.5, kappa_c = 0.5, kappa_s = 0.0)

First-order fully discrete scheme.

* `kind = :shift`: §3.3 (Theorem 3.5), `Q = η Δt ∂E φⁿ⁺¹`, η from (43), time step (44), (CFL).
* `kind = :dissipation`: §3.4 in the linearised variant (52), `Q = μ ∂E φⁿ⁺¹ - κ λ² ∂E ΔM φⁿ⁺¹`,
  `μ = η Δt + κ ρ̂ⁿ`, η from Corollary 3.17(i), time step of Remark 3.22, (CFLθ).
* `κⁿ = kappa_c Δx (max|uⁿ| + kappa_s)` (Remark 3.23; only for `:dissipation`).
"""
Base.@kwdef struct Scheme
    kind::Symbol = :shift
    theta::Float64 = 0.5
    kappa_c::Float64 = 0.5
    kappa_s::Float64 = 0.0
    function Scheme(kind, theta, kappa_c, kappa_s)
        kind in (:shift, :dissipation) ||
            throw(ArgumentError("unknown scheme kind :$kind; use :shift or :dissipation"))
        return new(kind, theta, kappa_c, kappa_s)
    end
end

"Local edge coefficient η_e: (43) for `:shift`, Corollary 3.17(i) for `:dissipation`."
function eta_coefficient(ρ, sch::Scheme; bc = :periodic)
    c = sch.kind == :shift ? 5 / 4 : 1 / (2 * (1 - sch.theta))
    return c .* edgemap(logmean, ρ; bc) .^ 2 ./ edgemap(harmonic, ρ; bc)
end

"Constant density-dissipation coefficient κⁿ (Remark 3.23); zero for `:shift`."
kappa_coefficient(u, dx, sch::Scheme) =
    sch.kind == :shift ? 0.0 : sch.kappa_c * dx * (maximum(abs, u) + sch.kappa_s)

# Sparse edge gradient D (rows: edges, columns: cells), so that divM = -Dᵀ and -ΔM = DᵀD.
# With walls the row of edge N is empty.
function gradient_matrix(N, dx; wall = false)
    E = wall ? N - 1 : N
    I = [1:E; 1:E]
    J = [1:E; [right(i, N) for i in 1:E]]
    V = [fill(-1 / dx, E); fill(1 / dx, E)]
    return sparse(I, J, V, N, N)
end

# Divergence of an :extrapolate edge array (N + 1 edges, the last one the left face of cell 1).
divergence_matrix(N, dx) = sparse([1:N; 1:N], [1:N; N + 1; 1:N-1], [fill(1 / dx, N); fill(-1 / dx, N)], N, N + 1)

"""
    solve_potential(φ, rhs, μ, κ, λ, dt, dx; bc = :periodic)

Newton's method for (52): `e^φ - λ² ΔM φ - Δt divM(μ ∂E φ) + κ λ² Δt ΔM² φ = rhs`, starting from `φ`.
The Jacobian is the SPD matrix of Remark 3.22, periodic corners included (none with walls); with
`:extrapolate`, μ has N + 1 edges and the matrix is not symmetric.
"""
function solve_potential(φ, rhs, μ, κ, λ, dt, dx; bc = :periodic, tol = 1e-12, maxit = 50)
    N = length(φ)
    D = gradient_matrix(N, dx; wall = bc == :wall)
    L = D' * D
    A = λ^2 * L + dt * (length(μ) == N ? D' * spdiagm(0 => μ) * D :
                        -divergence_matrix(N, dx) * spdiagm(0 => μ) * [D; D[N:N, :]])
    κ > 0 && (A += κ * λ^2 * dt * (L * L))
    φ = copy(φ)
    for _ in 1:maxit
        e = exp.(φ)
        R = e + A * φ - rhs
        J = A + spdiagm(0 => e)
        scale = maximum(abs, diag(J))          # stencil scale of Remark 4.13
        maximum(abs, R) <= 1e-2tol * scale * (1 + maximum(abs, φ)) && return φ
        δ = J \ R
        m = maximum(abs, δ)
        φ .-= δ ./ max(1, m)     # damp to |δφ| ≤ 1
        m <= tol * (1 + maximum(abs, φ)) && return φ
    end
    error("Newton iteration for the potential did not converge")
end

"Initial potential from ρ by (29)."
initial_potential(ρ, λ, dx; bc = :periodic) =
    solve_potential(log.(ρ), ρ, zeros(length(ρ)), 0.0, λ, 0.0, dx; bc)

"""
    step(ρ, u, φ, dt, dx, λ, η, κ; bc = :periodic)

One step of (27)–(29) with the given edge coefficients η and constant κ (κ = 0 gives §3.3);
with walls, F = G = 0 and ∂E φ = 0 on the walls; with `:extrapolate`, fluxes to ghost copies (grid.jl).
Returns `(ρ, u, φ, F)` at level n + 1, F being the mass flux used.
"""
function step(ρ, u, φ, dt, dx, λ, η, κ; bc = :periodic)
    N, wall = length(ρ), bc == :wall
    ρbar = edgemap(logmean, ρ; bc)                  # (17)
    conv = ρbar .* avg(u; bc)                       # ρ̄ {u}
    wall && (conv[end] = 0)
    μ = η .* dt
    κ > 0 && (μ = μ .+ κ .* widen(edgemap(logmean, exp.(φ)), bc))   # μ of (52), linearised ρ̂ⁿ
    φ1 = solve_potential(φ, ρ .- dt .* divm(conv, dx, N), μ, κ, λ, dt, dx; bc)
    dφ = dphi(φ1, dx, bc)
    Q = μ .* dφ
    κ > 0 && (Q = Q .- κ * λ^2 .* widen(grad(lap(φ1, dx; wall), dx; wall), bc))  # (46)
    F = conv .- Q                                   # (18)
    ρ1 = ρ .- dt .* divm(F, dx, N)                  # (27)
    uL, uR = edgestates(u; bc)
    G = ifelse.(F .>= 0, uL, uR) .* F               # (19)
    S = source(ρbar, dφ, N)                         # (15)
    m1 = ρ .* u .- dt .* divm(G, dx, N) .+ dt .* S  # (28)
    return ρ1, m1 ./ ρ1, φ1, F
end

"Source (15): `S_i = -(ρ̄ ∂E φ)_{i+1/2}/2 - (ρ̄ ∂E φ)_{i-1/2}/2` on `N` cells, from edge values ρ̄ and ∂E φ."
source(ρbar, dφ, N = length(ρbar)) = (a = ρbar .* dφ; [-(a[i] + lefte(a, i)) / 2 for i in 1:N])

# Outward flux of cell i through its right and left edges, convention (13).
outward(F, i) = (F[i], -lefte(F, i))

"""
    stable_dt(ρ, u, φ, dx, η, κ, sch; bc = :periodic)

Time step from (44) (`:shift`, with C_i = Δx Σ η|∂E φⁿ|) or Remark 3.22 times 0.9 (`:dissipation`);
no other bound (SPEC §5).
"""
function stable_dt(ρ, u, φ, dx, η, κ, sch::Scheme; bc = :periodic)
    N, wall = length(ρ), bc == :wall
    d = η .* abs.(dphi(φ, dx, bc))
    X = Inf
    if sch.kind == :shift
        a = abs.(edgemap(logmean, ρ; bc) .* avg(u; bc))
        wall && (a[end] = 0)
        for i in 1:N
            A = a[i] + lefte(a, i)
            C = dx * (d[i] + lefte(d, i))
            X = min(X, 0.4ρ[i] / (A + sqrt(A^2 + 0.8C * ρ[i])))
        end
        return X * dx
    end
    θ = sch.theta
    c = edgemap(logmean, ρ; bc) .* avg(u; bc) .- κ .* drho(ρ, dx, bc)
    wall && (c[end] = 0)
    for i in 1:N
        A = sum(max.(outward(c, i), 0))
        C = dx * (d[i] + lefte(d, i))
        X = min(X, 2θ * ρ[i] / (A + sqrt(A^2 + 4θ * C * ρ[i])))
    end
    return 0.9X * dx
end

"""
    cfl_ratio(ρ, F, dt, dx, sch)

Largest ratio of the left- to the right-hand side of (CFL) (`:shift`) or (CFLθ) (`:dissipation`),
evaluated with the computed flux F. The step is admissible iff the ratio is ≤ 1.
"""
function cfl_ratio(ρ, F, dt, dx, sch::Scheme)
    r = 0.0
    for i in eachindex(ρ)
        o = outward(F, i)
        if sch.kind == :shift
            r = max(r, dt / dx * sum(abs.(o)) / (ρ[i] / 5))
        else
            r = max(r, dt / dx * sum(max.(o, 0)) / (sch.theta * ρ[i]))
        end
    end
    return r
end
