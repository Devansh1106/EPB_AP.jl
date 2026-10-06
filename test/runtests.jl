using EPB_AP
using Test

@testset "EPB_AP.jl" begin
    include("test_grid.jl")
    include("test_scheme.jl")
    include("test_cases.jl")
    include("test_convergence.jl")
end
