"""
    Grid(a, b, N)

Uniform periodic grid of `N` cells on `[a, b]` with centres `x[i] = a + (i - 1/2) dx` (6).
"""
struct Grid
    a::Float64
    b::Float64
    N::Int
    dx::Float64
    x::Vector{Float64}
end

function Grid(a, b, N::Integer)
    dx = (b - a) / N
    return Grid(a, b, N, dx, [a + (i - 0.5) * dx for i in 1:N])
end

# Discrete operators (8). Edge e = i + 1/2 is stored at index i.

"Edge gradient `(∂E q)_{i+1/2} = (q_{i+1} - q_i)/dx`."
grad(q, dx) = [(q[right(i, length(q))] - q[i]) / dx for i in eachindex(q)]

"Cell divergence `(divM v)_i = (v_{i+1/2} - v_{i-1/2})/dx`."
divm(v, dx) = [(v[i] - v[left(i, length(v))]) / dx for i in eachindex(v)]

"Compact Laplacian `ΔM = divM ∘ ∂E`."
lap(q, dx) = divm(grad(q, dx), dx)

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
    energy(ρ, u, φ, λ, dx)

Discrete energy (14).
"""
function energy(ρ, u, φ, λ, dx)
    ekin = dx * sum(@. ρ * u^2 / 2)
    eel = dx * sum(@. exp(φ) * (φ - 1))
    efield = λ^2 / 2 * dx * sum(abs2, grad(φ, dx))
    return ekin + eel + efield
end
