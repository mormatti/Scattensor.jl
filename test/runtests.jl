using Test

# In package tests, the package is available on LOAD_PATH.
using Scattensor

@testset "Scattensor.jl" begin
    include("test_matrix_backend.jl")
    include("test_sectors.jl")
    include("test_localizability.jl")
    include("test_trap.jl")
    include("test_fits.jl")
    include("test_bulk_bands.jl")
    include("test_pair_bands.jl")
end