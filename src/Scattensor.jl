module Scattensor

    # Using  
    using ITensors, ITensorTDVP, ITensorGLMakie
    using LinearAlgebra, LinearSolve
    using Optim
    using Plots

    # Types
    include("types/exact_diag_system.jl")
    include("types/local_operator.jl")
    include("types/models.jl")

    # Algorithms
    include("algorithms/bloch.jl")
    include("algorithms/wannier.jl")
    include("algorithms/interpolation.jl")
    include("algorithms/translation_operator.jl")

    # Functions
    include("functions/notation.jl")
    include("functions/operators.jl")
    include("functions/entanglement_entropy.jl")
    include("functions/tn_evolution.jl")

    # Test
    include("test.jl")

end # module Scattensor