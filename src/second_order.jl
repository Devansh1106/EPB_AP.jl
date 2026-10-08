"""
    SecondOrder(; limiter = :minmod, c_D = 1.0, cfl = 0.4)

Second-order scheme of §4 (SPEC §10): MUSCL reconstruction (54)–(55) on (ρ, v) with
`limiter = :minmod` or `:none` (central slopes), fluxes (56)–(63) with the upwind momentum flux (61),
density dissipation `D = c_D Δx² ∂E ρ` (72), the ARS(2,2,2) stage sequence (76)–(87) and the time step
min((44), (91)) with `cfl` in (91).
"""
Base.@kwdef struct SecondOrder
    limiter::Symbol = :minmod
    c_D::Float64 = 1.0
    cfl::Float64 = 0.4
    function SecondOrder(limiter, c_D, cfl)
        limiter in (:minmod, :none) ||
            throw(ArgumentError("unknown limiter :$limiter; use :minmod or :none"))
        return new(limiter, c_D, cfl)
    end
end

"Minmod of two slopes (§4.1)."
minmod(a, b) = a * b <= 0 ? zero(a) : (abs(a) < abs(b) ? a : b)

"Cell slopes σ of (54); zero in cells 1 and N with walls."
function slopes(w, dx, limiter; wall = false)
    N = length(w)
    σ = map(1:N) do i
        l, r = w[i] - w[left(i, N)], w[right(i, N)] - w[i]
        limiter == :minmod ? minmod(l, r) / dx : (l + r) / 2dx
    end
    wall && (σ[1] = σ[N] = 0)
    return σ
end

"Edge states `(w⁻, w⁺)` of (55) at e = i + 1/2 (index i), and the slopes."
function reconstruct(w, dx, limiter; wall = false)
    σ = slopes(w, dx, limiter; wall)
    N = length(w)
    wm = [w[i] + dx / 2 * σ[i] for i in 1:N]
    wp = [w[right(i, N)] - dx / 2 * σ[right(i, N)] for i in 1:N]
    return wm, wp, σ
end

"Thrown when a stage or reconstructed density is not positive."
struct NonPositiveDensity <: Exception end

"""
    faces(ρ, q, dx, sch; wall = false)

Edge quantities of the state (ρ, q): interface density ρ̄ = LM(ρ⁻, ρ⁺) (56), explicit flux
E = ρ̄ {v} (58), dissipation D (72), reconstructed velocities v∓ and the velocity slopes.
"""
function faces(ρ, q, dx, sch::SecondOrder; wall = false)
    minimum(ρ) > 0 || throw(NonPositiveDensity())
    v = q ./ ρ
    ρm, ρp, _ = reconstruct(ρ, dx, sch.limiter; wall)
    min(minimum(ρm), minimum(ρp)) > 0 || throw(NonPositiveDensity())
    vm, vp, σv = reconstruct(v, dx, sch.limiter; wall)
    ρbar = logmean.(ρm, ρp)
    E = ρbar .* avg(v)
    D = sch.c_D * dx^2 .* grad(ρ, dx; wall)
    wall && (E[end] = 0)
    return (ρbar = ρbar, E = E, D = D, vm = vm, vp = vp, σv = σv, v = v)
end

"Upwind momentum flux (61): `G = v⁻ F⁺ + v⁺ F⁻`."
upwind_flux(F, f) = f.vm .* max.(F, 0) .+ f.vp .* min.(F, 0)

"""
    semidiscrete_rhs(ρ, q, φ, dx, sch; wall = false) -> (ρ̇, q̇, F, f)

Right-hand side of (65)–(66) with `F = E - D` (no semi-implicit part); `f` are the [`faces`](@ref).
"""
function semidiscrete_rhs(ρ, q, φ, dx, sch::SecondOrder; wall = false)
    f = faces(ρ, q, dx, sch; wall)
    F = f.E .- f.D
    S = source(f.ρbar, grad(φ, dx; wall))
    return -divm(F, dx), -divm(upwind_flux(F, f), dx) .+ S, F, f
end

# ARS(2,2,2) pair (73)–(74).
const GAMMA = 1 - 1 / sqrt(2)
const DELTA = 1 - 1 / (2GAMMA)

"Stage mass flux (60) with the semi-implicit part (59): `F = E - γΔt ρ̄ ∂E φ - D`."
stage_flux(f, φ, dt, dx; wall = false) = f.E .- GAMMA * dt .* f.ρbar .* grad(φ, dx; wall) .- f.D

"Elliptic solve (88), `e^φ - λ² ΔM φ - γ²Δt² divM(ρ̄ ∂E φ) = ρ̂`, by the Newton solver of §6 from `φ`."
stage_potential(φ, ρhat, f, dt, λ, dx; wall = false) =
    solve_potential(φ, ρhat, GAMMA^2 * dt .* f.ρbar, 0.0, λ, dt, dx; wall)

"""
    step(ρ, u, φ, dt, dx, λ, sch::SecondOrder; wall = false)

One step of the stage sequence (76)–(87) (SPEC §10). Returns `(ρ, u, φ)` at level n + 1; throws
`NonPositiveDensity` if a stage or reconstructed density is not positive.
"""
function step(ρ, u, φ, dt, dx, λ, sch::SecondOrder; wall = false)
    γ, δ = GAMMA, DELTA
    q = ρ .* u
    div(F) = divm(F, dx)
    # Stage 1 at (ρⁿ, qⁿ, φⁿ).
    f1 = faces(ρ, q, dx, sch; wall)
    F1 = f1.E .- f1.D
    C1 = div(upwind_flux(F1, f1))
    S1 = source(f1.ρbar, grad(φ, dx; wall))
    # Stage 2.
    ρE2 = ρ .- γ * dt .* div(F1)                                        # (76)
    qt2 = q .- γ * dt .* C1                                             # (77) force-free
    qE2 = qt2 .+ γ * dt .* S1                                           # (78) force-carrying
    f2 = faces(ρE2, qt2, dx, sch; wall)
    φ2 = stage_potential(φ, ρ .- γ * dt .* div(f2.E .- f2.D), f2, dt, λ, dx; wall)   # (79)
    F2 = stage_flux(f2, φ2, dt, dx; wall)                               # (80)
    S2 = source(f2.ρbar, grad(φ2, dx; wall))
    # Stage 3.
    C2 = div(upwind_flux(F2, faces(ρE2, qE2, dx, sch; wall)))           # (81)
    ρE3 = ρ .- dt .* div(δ .* F1 .+ (1 - δ) .* F2)                      # (82)
    qh3 = q .- dt .* (δ .* C1 .+ (1 - δ) .* C2) .+ (1 - γ) * dt .* S2    # (83)
    f3 = faces(ρE3, qh3, dx, sch; wall)                                 # (84)
    φ3 = stage_potential(φ2, ρ .- dt .* div((1 - γ) .* F2 .+ γ .* (f3.E .- f3.D)), f3, dt, λ, dx; wall)  # (85)
    F3 = stage_flux(f3, φ3, dt, dx; wall)
    ρ1 = ρ .- dt .* div((1 - γ) .* F2 .+ γ .* F3)                       # (86)
    minimum(ρ1) > 0 || throw(NonPositiveDensity())
    q1 = qh3 .+ γ * dt .* source(f3.ρbar, grad(φ3, dx; wall))           # (87)
    return ρ1, q1 ./ ρ1, φ3
end

"""
    stable_dt(ρ, u, φ, dx, sch::SecondOrder; wall = false)

Time step of §4.11: the minimum of (44) and (91), both with η of (43), `cfl = sch.cfl` in (91).
"""
function stable_dt(ρ, u, φ, dx, sch::SecondOrder; wall = false)
    η = eta_coefficient(ρ, Scheme())
    dt44 = stable_dt(ρ, u, φ, dx, η, 0.0, Scheme(); wall)
    wall && (η[end] = 0)
    jump = abs.(grad(φ, dx; wall))                      # |φ_{i+1} - φ_i| / Δx
    N = length(ρ)
    s = maximum(abs(u[i]) + max(η[i], η[left(i, N)]) * max(jump[i], jump[left(i, N)]) for i in 1:N)
    return s > 0 ? min(dt44, sch.cfl * dx / s) : dt44
end
