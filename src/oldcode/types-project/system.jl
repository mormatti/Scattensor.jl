"Struct for a quantum many body system on a uniform, on a finite one-dimensional 
lattice with periodic boundary conditions. This system uses matrices, hence
it is not suitable for large systems. It is implemented to perform exact
diagonalization algorithms."

# TYPE DEFINITION

mutable struct ExactDiagSystem
    # L, the number of sites of the chain
    size                        :: Int64
    # d, The local dimension
    local_dimension             :: Int64
    # 𝐓, the translation operator
    translation_operator        :: ExactDiagMatrixFormOperator
    # 𝐇, the Hamiltonian operator
    hamiltonian_operator        :: ExactDiagMatrixFormOperator
    # 𝛙[i], the eigenstates of the system
    eigenstates                 :: Union{ExactDiagStateSet, Missing}
    # ℰ₀, the groundstate of the system
    groundstate                 :: Union{ExactDiagState, Missing}
    # The Bloch states of the first band
    first_band_states           :: Union{Vector{ExactDiagState}, Missing}
end
# export ExactDiagSystem later

function ExactDiagSystem(size::Int64, local_dimension::Int64, hamiltonian_operator::Union{ExactDiagOperator, Missing})
    return ExactDiagSystem(size, local_dimension, translation_operator, missing, missing, missing, missing)
end

size(𝒮::ExactDiagSystem) = 𝒮.size

local_dimension(𝒮::ExactDiagSystem) = 𝒮.local_dimension

translation_operator(𝒮::ExactDiagSystem) = 𝒮.translation_operator

hamiltonian_operator(𝒮::ExactDiagSystem) = 𝒮.hamiltonian_operator

eigenstates(𝒮::ExactDiagSystem) = 𝒮.eigenstates

groundstate(𝒮::ExactDiagSystem) = 𝒮.groundstate

first_band_states(𝒮::ExactDiagSystem) = 𝒮.first_band_states


# CONSTRUCTORS

function ExactDiagSystem(L::Integer, d::Integer)::ExactDiagSystem
    ∅ = missing
    𝐓::Matrix{ComplexF64} = translation_operator(L, d)
    𝒮 = ExactDiagSystem(L, d, 𝐓, ∅, ∅, ∅, ∅)
    return 𝒮
end
export ExactDiagSystem


# STRUCT FUNCTIONS

presence_𝐇(𝒮::ExactDiagSystem) = !ismissing(𝒮.hamiltonian_operator)
export presence_𝐇

presence_eigenstates(𝒮::ExactDiagSystem) = !ismissing(𝒮.eigenstates)
export presence_eigenstates

function check_translational_invariance(𝒮::ExactDiagSystem)
    @assert presence_𝐇(𝒮) "𝐇 not already inserted."
    𝐇::Matrix{ComplexF64} = 𝒮.hamiltonian_operator
    𝐓::Matrix{ComplexF64} = 𝒮.translation_operator
    if (𝐇 * 𝐓 == 𝐓 * 𝐇)
        @warn "Found hamiltonian not translational invariant."
        𝒮.translational_invariance = false
    else
        @debug "Translational invariance of 𝐇 checked."
        𝒮.translational_invariance = true
    end
end
export check_translational_invariance