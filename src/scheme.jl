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
end

"Local edge coefficient η_e: (43) for `:shift`, Corollary 3.17(i) for `:dissipation`."
function eta_coefficient(ρ, sch::Scheme)
    c = sch.kind == :shift ? 5 / 4 : 1 / (2 * (1 - sch.theta))
    return c .* edgemap(logmean, ρ) .^ 2 ./ edgemap(harmonic, ρ)
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

"""
    solve_potential(φ, rhs, μ, κ, λ, dt, dx; wall = false)

Newton's method for (52): `e^φ - λ² ΔM φ - Δt divM(μ ∂E φ) + κ λ² Δt ΔM² φ = rhs`, starting from `φ`.
The Jacobian is the SPD matrix of Remark 3.22, periodic corners included (none with walls).
"""
function solve_potential(φ, rhs, μ, κ, λ, dt, dx; wall = false, tol = 1e-12, maxit = 50)
    N = length(φ)
    D = gradient_matrix(N, dx; wall)
    L = D' * D
    A = λ^2 * L + dt * (D' * spdiagm(0 => μ) * D)
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
initial_potential(ρ, λ, dx; wall = false) =
    solve_potential(log.(ρ), ρ, zeros(length(ρ)), 0.0, λ, 0.0, dx; wall)

"""
    step(ρ, u, φ, dt, dx, λ, η, κ; wall = false)

One step of (27)–(29) with the given edge coefficients η and constant κ (κ = 0 gives §3.3);
with `wall`, F = G = 0 and ∂E φ = 0 on the walls.
Returns `(ρ, u, φ, F)` at level n + 1, F being the mass flux used.
"""
function step(ρ, u, φ, dt, dx, λ, η, κ; wall = false)
    N = length(ρ)
    ρbar = edgemap(logmean, ρ)                      # (17)
    conv = ρbar .* avg(u)                           # ρ̄ {u}
    wall && (conv[end] = 0)
    μ = η .* dt
    κ > 0 && (μ = μ .+ κ .* edgemap(logmean, exp.(φ)))   # μ of (52), linearised ρ̂ⁿ
    φ1 = solve_potential(φ, ρ .- dt .* divm(conv, dx), μ, κ, λ, dt, dx; wall)
    dφ = grad(φ1, dx; wall)
    Q = μ .* dφ
    κ > 0 && (Q = Q .- κ * λ^2 .* grad(lap(φ1, dx; wall), dx; wall))  # (46)
    F = conv .- Q                                   # (18)
    ρ1 = ρ .- dt .* divm(F, dx)                     # (27)
    G = [F[e] >= 0 ? u[e] * F[e] : u[right(e, N)] * F[e] for e in 1:N]   # (19)
    a = ρbar .* dφ
    S = [-(a[i] + a[left(i, N)]) / 2 for i in 1:N]  # (15)
    m1 = ρ .* u .- dt .* divm(G, dx) .+ dt .* S     # (28)
    return ρ1, m1 ./ ρ1, φ1, F
end

# Outward flux of cell i through edge e = i (right) and e = i-1 (left), convention (13).
outward(F, i) = (F[i], -F[left(i, length(F))])

"""
    stable_dt(ρ, u, φ, dx, η, κ, sch; wall = false)

Time step from (44) (`:shift`, with C_i = Δx Σ η|∂E φⁿ|) or Remark 3.22 times 0.9 (`:dissipation`);
no other bound (SPEC §5).
"""
function stable_dt(ρ, u, φ, dx, η, κ, sch::Scheme; wall = false)
    N = length(ρ)
    d = η .* abs.(grad(φ, dx; wall))
    X = Inf
    if sch.kind == :shift
        a = abs.(edgemap(logmean, ρ) .* avg(u))
        wall && (a[end] = 0)
        for i in 1:N
            A = a[i] + a[left(i, N)]
            C = dx * (d[i] + d[left(i, N)])
            X = min(X, 0.4ρ[i] / (A + sqrt(A^2 + 0.8C * ρ[i])))
        end
        return X * dx
    end
    θ = sch.theta
    c = edgemap(logmean, ρ) .* avg(u) .- κ .* grad(ρ, dx; wall)
    wall && (c[end] = 0)
    for i in 1:N
        A = sum(max.(outward(c, i), 0))
        C = dx * (d[i] + d[left(i, N)])
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
