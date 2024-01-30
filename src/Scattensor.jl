module Scattensor

    # Using  
    using ITensors, ITensorTDVP, ITensorGLMakie
    using LinearAlgebra, LinearSolve
    using SparseArrays
    using Optim
    using Plots

    # Types
    include("types/exact_diag/system.jl")
    include("types/exact_diag/operator.jl")
    include("types/exact_diag/local_operator.jl")
    include("types/exact_diag/state.jl")
    include("types/exact_diag/states_set.jl")
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