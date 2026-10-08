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
