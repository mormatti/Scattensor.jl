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

    # Abstract functions
    include("abstract_functions/operator_identity.jl")
    include("abstract_functions/kron_power.jl")

    # Utils
    include("utils/periodic_modulus.jl")
    include("utils/plot_complex_colormap.jl")
    include("utils/print_with_color.jl")
    include("utils/rgb_from_angle.jl")
    include("utils/rgb_from_complex.jl")
    include("utils/rgb_from_hex.jl")

    # Matrix extensions
    include("matrix_extensions/hilbspace_warning.jl")
    include("matrix_extensions/kron_power.jl")
    include("matrix_extensions/mathematica_format.jl")
    include("matrix_extensions/operator_identity.jl")
    include("matrix_extensions/operator_reflection.jl")
    include("matrix_extensions/operator_translation.jl")
    include("matrix_extensions/product_locals.jl")
    include("matrix_extensions/summation_local.jl")

    # ITensor extensions
    include("itensor_extensions/apply_reflection.jl")
    include("itensor_extensions/apply_translation.jl")
    include("itensor_extensions/entanglement_entropy.jl")
    include("itensor_extensions/insert_local.jl")
    include("itensor_extensions/local_expvals.jl")
    include("itensor_extensions/mpo_from_matrix.jl")
    include("itensor_extensions/mps_from_vector.jl")
    include("itensor_extensions/operator_identity.jl")
    include("itensor_extensions/operator_translation.jl")
    include("itensor_extensions/partial_trace.jl")
    include("itensor_extensions/kron.jl")
    include("itensor_extensions/product_inner.jl")
    include("itensor_extensions/product_matricial.jl")
    include("itensor_extensions/product_outer.jl")
    include("itensor_extensions/substitute_siteinds.jl")
    include("itensor_extensions/summation_local.jl")
    include("itensor_extensions/tdvp_time_evolution.jl")

    # Scattensor functions
    include("scattensor_functions/blochstate.jl")
    include("scattensor_functions/compute_smatrix.jl")
    include("scattensor_functions/diagonalization_HU.jl")
    include("scattensor_functions/dispersion_relation.jl")
    include("scattensor_functions/fourier_transform.jl")
    include("scattensor_functions/get_firstband.jl")
    include("scattensor_functions/get_groundstate.jl")
    include("scattensor_functions/get_statesabove.jl")
    include("scattensor_functions/get_statesbelow.jl")
    include("scattensor_functions/local_exp_value.jl")
    include("scattensor_functions/plot_disprel.jl")
    include("scattensor_functions/pop_groundstate.jl")
    include("scattensor_functions/wannier_symmetric.jl")

end # module Scattensor