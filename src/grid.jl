"""
    Grid(a, b, N; bc = :periodic)

Uniform grid of `N` cells on `[a, b]` with centres `x[i] = a + (i - 1/2) dx` (6).
`bc = :periodic` or `:wall` (walls at a and b, SPEC §7).
"""
struct Grid
    a::Float64
    b::Float64
    N::Int
    dx::Float64
    x::Vector{Float64}
    bc::Symbol
end

function Grid(a, b, N::Integer; bc = :periodic)
    bc in (:periodic, :wall) || throw(ArgumentError("bc must be :periodic or :wall"))
    dx = (b - a) / N
    return Grid(a, b, N, dx, [a + (i - 0.5) * dx for i in 1:N], bc)
end

# Discrete operators (8). Edge e = i + 1/2 is stored at index i. With walls, the edge stored at
# index N stands for both walls (1/2 and N + 1/2): every edge quantity there is zero (no flux,
# ∂E φ = 0), so divM at cells 1 and N sees the wall.

"Edge gradient `(∂E q)_{i+1/2} = (q_{i+1} - q_i)/dx`; zero at the walls if `wall`."
function grad(q, dx; wall = false)
    g = [(q[right(i, length(q))] - q[i]) / dx for i in eachindex(q)]
    wall && (g[end] = 0)
    return g
end

"Cell divergence `(divM v)_i = (v_{i+1/2} - v_{i-1/2})/dx`."
divm(v, dx) = [(v[i] - v[left(i, length(v))]) / dx for i in eachindex(v)]

"Compact Laplacian `ΔM = divM ∘ ∂E` (homogeneous Neumann at the walls if `wall`)."
lap(q, dx; wall = false) = divm(grad(q, dx; wall), dx)

"Edge function from a cell function via `f(q_i, q_{i+1})`."
edgemap(f, q) = [f(q[i], q[right(i, length(q))]) for i in eachindex(q)]

"Arithmetic edge average `{u}_e` (12)."
avg(u) = edgemap((a, b) -> (a + b) / 2, u)

"""
    logmean(a, b)

Logarithmic mean (9), evaluated with a series near `a = b` to avoid cancellation.
"""
function logmean(a, b)
    f = (a - b) / (a + b)
    v = f^2
    if v < 1e-2
        F = 1 + v * (1 / 3 + v * (1 / 5 + v * (1 / 7 + v * (1 / 9 + v * (1 / 11 + v / 13)))))
        return (a + b) / (2F)
    end
    return (b - a) / (log(b) - log(a))
end

"Harmonic edge mean (11)."
harmonic(a, b) = 2a * b / (a + b)

"""
    energy(ρ, u, φ, λ, dx; wall = false)

Discrete energy (14).
"""
function energy(ρ, u, φ, λ, dx; wall = false)
    ekin = dx * sum(@. ρ * u^2 / 2)
    eel = dx * sum(@. exp(φ) * (φ - 1))
    efield = λ^2 / 2 * dx * sum(abs2, grad(φ, dx; wall))
    return ekin + eel + efield
end
