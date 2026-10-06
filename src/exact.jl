# Exact and reference solutions (SPEC §7). Each returns (ρ, u, φ) at the points x and time t.

"""
    SolitonProfile(c, ε, ξmax; h = 1e-3)

Ion-acoustic travelling wave of speed c (§5.1): ε² φ'' = e^φ - ρ(φ), ρ = c/√(c² - 2φ),
integrated by RK4 from the crest (φ_max, 0), where V(φ_max) = 0.
"""
struct SolitonProfile
    c::Float64
    h::Float64
    phi::Vector{Float64}
    dphi::Vector{Float64}
end

function SolitonProfile(c, ε, ξmax; h = 1e-3)
    V(φ) = exp(φ) - 1 + c * (sqrt(c^2 - 2φ) - c)
    lo, hi = 1e-3, c^2 / 2                  # V(lo) > 0 > V(hi) for c > 1
    for _ in 1:200
        mid = (lo + hi) / 2
        V(mid) > 0 ? (lo = mid) : (hi = mid)
    end
    f(φ, ψ) = (ψ, (exp(φ) - c / sqrt(c^2 - 2φ)) / ε^2)
    n = ceil(Int, ξmax / h) + 1
    φ, ψ = zeros(n + 1), zeros(n + 1)
    φ[1] = lo
    for k in 1:n
        k1 = f(φ[k], ψ[k])
        k2 = f(φ[k] + h / 2 * k1[1], ψ[k] + h / 2 * k1[2])
        k3 = f(φ[k] + h / 2 * k2[1], ψ[k] + h / 2 * k2[2])
        k4 = f(φ[k] + h * k3[1], ψ[k] + h * k3[2])
        φ[k+1] = φ[k] + h / 6 * (k1[1] + 2k2[1] + 2k3[1] + k4[1])
        ψ[k+1] = ψ[k] + h / 6 * (k1[2] + 2k2[2] + 2k3[2] + k4[2])
    end
    return SolitonProfile(c, h, φ, ψ)
end

"Potential of the profile at distance ξ from the crest (cubic Hermite interpolation)."
function (p::SolitonProfile)(ξ)
    s = abs(ξ) / p.h
    k = min(floor(Int, s), length(p.phi) - 2)
    τ = s - k
    y0, y1, d0, d1 = p.phi[k+1], p.phi[k+2], p.h * p.dphi[k+1], p.h * p.dphi[k+2]
    return (2τ^3 - 3τ^2 + 1) * y0 + (τ^3 - 2τ^2 + τ) * d0 + (-2τ^3 + 3τ^2) * y1 + (τ^3 - τ^2) * d1
end

"Exact soliton on the periodic interval [a, b], crest at the centre at t = 0."
function soliton_exact(p::SolitonProfile, a, b)
    L, x0, c = b - a, (a + b) / 2, p.c
    return function (x, t)
        φ = [p(mod(xi - x0 - c * t + L / 2, L) - L / 2) for xi in x]
        s = sqrt.(c^2 .- 2φ)
        return c ./ s, c .- s, φ
    end
end

"Isothermal simple wave u - log ρ = 0 with u⁰ = A sin 2πx, constant along ẋ = u + 1 (§5.1)."
function simplewave_exact(x, t; A = 0.2)
    u = map(x) do xi
        ξ = xi - t
        for _ in 1:100
            g = ξ + (A * sin(2π * ξ) + 1) * t - xi
            ξ -= g / (1 + 2π * A * t * cos(2π * ξ))
            abs(g) < 1e-15 && break
        end
        A * sin(2π * ξ)
    end
    return exp.(u), u, u
end

"Middle velocity u_m and shock speed u_s of the ICE Riemann problem, [1, (7.9)–(7.10)]."
function riemann_speeds(nr)
    g(um) = (1 - nr * exp(um)) * (um^2 - 2um - 2log(nr)) - 2um^2
    lo, hi = 0.0, -log(nr)                   # g(lo) > 0 > g(hi)
    for _ in 1:200
        mid = (lo + hi) / 2
        g(mid) > 0 ? (lo = mid) : (hi = mid)
    end
    um = (lo + hi) / 2
    return um, um / (1 - nr * exp(um))
end

"Self-similar ICE solution [1, (7.8)]; φ = log ρ (quasineutral limit)."
function riemann_exact(x, t; nr = 0.5)
    um, us = riemann_speeds(nr)
    ρu = map(x) do xi
        xi <= -t ? (1.0, 0.0) :
        xi <= (um - 1) * t ? (exp(-xi / t - 1), xi / t + 1) :
        xi <= us * t ? (exp(-um), um) : (nr, 0.0)
    end
    ρ = first.(ρu)
    return ρ, last.(ρu), log.(ρ)
end

"Shock position by equal areas over the cells with centres in [½(u_m - 1 + u_s)t, xb] (§5.5)."
function shock_position(x, ρ, t, dx; nr = 0.5, xb = 60.0)
    um, us = riemann_speeds(nr)
    xa = (um - 1 + us) * t / 2
    idx = findall(xi -> xa <= xi <= xb, x)
    return x[first(idx)] - dx / 2 + dx * sum(ρ[idx] .- nr) / (exp(-um) - nr)
end
