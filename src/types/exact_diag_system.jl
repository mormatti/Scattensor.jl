"Struct for a quantum many body system on a uniform, on a finite one-dimensional 
lattice with periodic boundary conditions. This system uses matrices, hence
it is not suitable for large systems. It is implemented to perform exact
diagonalization algorithms."

# TYPE DEFINITION

mutable struct ExactDiagSystem
    # L, the number of sites of the chain
    system_size                 :: Int64
    # d, The local dimension
    local_dimension             :: Int64
    # 𝐇, the Hamiltonian operator
    hamiltonian_operator        :: Union{Matrix{ComplexF64}, Missing}
    # 𝐓, the translation operator
    translation_operator        :: Matrix{ComplexF64}
    # Whether [H,T] = 0
    translational_invariance    :: Union{Bool, Missing} 
    # 𝛙[i], the eigenstates of the system
    eigenstates                 :: Union{Vector{Vector{ComplexF64}}, Missing}
    # ℰ[i], the eigenenergies of the system
    eigenstates_energies        :: Union{Vector{Float64}, Missing}
    # k[i], the eigenmomenta of the system
    eigenstates_momenta         :: Union{Vector{Float64}, Missing}
    # ℰ₀, the groundstate of the system
    groundstate                 :: Union{Vector{ComplexF64}, Missing}
    # The groundstate energy of the system
    groundstate_energy          :: Union{Float64, Missing}
    # The Bloch states of the first band
    first_band_states           :: Union{Vector{Vector{ComplexF64}}, Missing}
end
# export later


# CONSTRUCTORS

function ExactDiagSystem(L::Integer, d::Integer)::ExactDiagSystem
    𝐓::Matrix{ComplexF64} = translation_operator(L, d)
    𝒮 = ExactDiagSystem(L,d,missing,𝐓,missing,missing,missing,missing,missing,missing,missing)
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