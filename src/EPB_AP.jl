module EPB_AP

using LinearAlgebra
using SparseArrays

include("boundary.jl")
include("grid.jl")
include("scheme.jl")

export Grid, energy, Scheme

end
