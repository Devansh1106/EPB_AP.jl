"""
    Grid(a, b, N; bc = :periodic)

Uniform grid of `N` cells on `[a, b]` with centres `x[i] = a + (i - 1/2) dx` (6).
`bc = :periodic`, `:wall` (walls at a and b), `:extrapolate` (ρ, u extrapolated, φ periodic) or `:neumann`
(ρ extrapolated, u = 0 in the left ghost and extrapolated at the right, φ homogeneous Neumann), SPEC §7.
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
    bc in (:periodic, :wall, :extrapolate, :neumann) ||
        throw(ArgumentError("bc must be :periodic, :wall, :extrapolate or :neumann"))
    dx = (b - a) / N
    return Grid(a, b, N, dx, [a + (i - 0.5) * dx for i in 1:N], bc)
end

# Discrete operators (8). Edge e = i + 1/2 is stored at index i. With walls, the edge stored at
# index N stands for both walls (1/2 and N + 1/2): every edge quantity there is zero (no flux,
# ∂E φ = 0), so divM at cells 1 and N sees the wall. With :extrapolate, edge arrays of the fluxes have
# N + 1 entries: edge N is the right face of cell N and edge N + 1 the left face of cell 1, each between
# the cell and its ghost copy, while φ stays periodic (its edge gradient at N + 1 repeats that at N).
# `:neumann` has the layout of :extrapolate, a zero velocity in the left ghost and ∂E φ = 0 at both
# boundary faces, as with walls. The left edge of cell 1 is always the last entry.

"Flux edge arrays have N + 1 entries, the last two edges facing ghost cells."
ghosted(bc) = bc in (:extrapolate, :neumann)

"φ satisfies homogeneous Neumann conditions at both ends."
neumann_phi(bc) = bc in (:wall, :neumann)

"Velocity in the left ghost cell (`nothing`: a copy of cell 1)."
vghost(bc) = bc == :neumann ? 0.0 : nothing

"Value of the edge array `v` at the left edge of cell `i`."
@inline lefte(v, i) = v[i == 1 ? lastindex(v) : i - 1]

"Edge gradient `(∂E q)_{i+1/2} = (q_{i+1} - q_i)/dx`; zero at the walls if `wall`."
function grad(q, dx; wall = false)
    g = [(q[right(i, length(q))] - q[i]) / dx for i in eachindex(q)]
    wall && (g[end] = 0)
    return g
end

"Cell divergence `(divM v)_i = (v_{i+1/2} - v_{i-1/2})/dx` on `N` cells."
divm(v, dx, N = length(v)) = [(v[i] - lefte(v, i)) / dx for i in 1:N]

"Compact Laplacian `ΔM = divM ∘ ∂E` (homogeneous Neumann at the walls if `wall`)."
lap(q, dx; wall = false) = divm(grad(q, dx; wall), dx)

"""
    edgestates(r, l = r; bc = :periodic, ghost = nothing)

States on the two sides of every edge, from the values `r` (`l`) that each cell gives its right (left)
face; with `:extrapolate` or `:neumann` the ghost beyond a boundary face repeats the cell inside, except
that the left ghost has the value `ghost` if given.
"""
function edgestates(r, l = r; bc = :periodic, ghost = nothing)
    N = length(r)
    ghosted(bc) || return r, l[right.(1:N, N)]
    g = ghost === nothing ? l[1] : ghost
    return [r; g], [l[2:N]; r[N]; l[1]]
end

"Edge function from a cell function via `f(q_i, q_{i+1})`."
edgemap(f, q; bc = :periodic, ghost = nothing) = f.(edgestates(q; bc, ghost)...)

"Arithmetic edge average `{u}_e` (12)."
avg(u; bc = :periodic, ghost = nothing) = edgemap((a, b) -> (a + b) / 2, u; bc, ghost)

"Edge array of a periodic or Neumann quantity in the layout of `bc` (edge N repeated with ghosts)."
widen(g, bc) = ghosted(bc) ? [g; g[end]] : g

"Edge gradient of φ for the fluxes: periodic with `:extrapolate`, Neumann with walls and `:neumann`."
dphi(φ, dx, bc) = widen(grad(φ, dx; wall = neumann_phi(bc)), bc)

"Edge gradient of a hyperbolic variable: through the ghosts with `:extrapolate`, `:neumann` (zero at the boundary faces)."
drho(ρ, dx, bc) = ghosted(bc) ? edgemap((a, b) -> (b - a) / dx, ρ; bc) : grad(ρ, dx; wall = bc == :wall)

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
