module EPB_AP

using LinearAlgebra
using SparseArrays
using Printf

include("boundary.jl")
include("grid.jl")
include("scheme.jl")
include("solve.jl")
include("testcases.jl")
include("output.jl")

export Grid, energy, Scheme, TimeStep, solve, testcase, run_case, read_csv, write_csv

end
