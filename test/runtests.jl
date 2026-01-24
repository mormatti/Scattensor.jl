using Test

# In package tests, the package is available on LOAD_PATH.
using Scattensor

@testset "Scattensor.jl" begin
    include("test_matrix_backend.jl")
end