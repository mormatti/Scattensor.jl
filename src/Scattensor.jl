module Scattensor

    # Usings 
    using ITensors, ITensorMPS
    using LinearAlgebra
    using SparseArrays
    using Optim
    using Plots
    using PlotlyJS
    using KrylovKit
    using Colors
    using Logging
    using LaTeXStrings

    # Utils
    include("utils/notations.jl")
    include("utils/outputs.jl")

    # Types
    include("types/localoperator.jl")
    include("types/blochstate.jl")

    # Conversions
    include("utils/conversions.jl")

    # Operators
    include("operators/identity.jl")
    include("operators/translation.jl")
    include("operators/reflection.jl")

    # Functions
    include("functions/localops.jl")
    include("functions/products.jl")
    include("functions/diagonalization.jl")
    include("functions/entropy.jl")
    include("functions/partialtrace.jl")
    include("functions/evolution.jl")

    # Scattering
    include("scattering/disprel.jl")
    include("scattering/wannier.jl")
    include("functions/smatrix.jl")

end # module Scattensor