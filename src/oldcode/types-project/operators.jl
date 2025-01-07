"""
Generates the whole matrix corresponding to the action of a product of local operators.

Inputs:
- `𝒮` is the quantum system;
- args is a list of pairs (𝐚, j) where 𝐚 is the local operator written in the local 
space (small matrix) and j is the position of the local operator 𝐚.
"""

function product_local_operators(𝒮::ExactDiagSystem, args::Vararg{LocalOperator})::Matrix{ComplexF64}
    L = 𝒮.system_size
    d = 𝒮.local_dimension
    N = d^L
    # We assert that all the matrices have the same dimension of the local space
    for arg in args
        @assert (size(arg.matrix) == (d, d)) "The local operator must have the same dimension of the local space."
    end
    # For each arg.second we take the mod L using ↻L
    for arg in args
        arg.position = arg.position ↻ L
    end
    # We consider the dxd identity matrix
    𝟙 = Matrix{ComplexF64}(I, d, d)
    # We create a Vector 𝒱 of Vector of Matrix{ComplexF64, where 𝒱[i] is the vector 
    # which contains all the operators of the i-th site
    𝒱 = []
    for i in 1:L
        𝒱i = []
        for arg in args
            if arg[2] == i
                push!(𝒱i, arg[1])
            end
        end
        # If 𝒱i has no elements, we add the identity
        if length(𝒱i) == 0
            push!(𝒱i, 𝟙)
        end
        push!(𝒱, 𝒱i)
    end
    # We substitute each 𝒱[i] with a Matrix{ComplexF64} which corresponds to the product
    # of all the operators in 𝒱[i]
    for i in 1:L
        𝒱[i] = reduce(*, 𝒱[i])
    end
    @assert length(𝒱) == L "Something went wrong in function local_operator."
    # We generate a tensor product of L dxd matrices: each element of this tensor from 𝒱[i]
    # is a factor of the tensor product (kronecker product ⊗)
    𝐀 = reduce(⊗, 𝒱)
    @assert size(𝐀) == (N, N) "Something went wrong in function local_operator."
    # We return the matrix 𝐀
    return 𝐀
end
export product_local_operators


"""
Generates the Hamiltonian matrix of the Ising model with transverse and longitudinal fields.

Inputs:
- `J` is the coupling constant of the spins;
- `hˣ` is the transverse field;
- `hᶻ` is the longitudinal field;
- `L` is the number of sites of the chain.
"""
function ising_hamiltonian(
    𝒮::ExactDiagSystem, # The quantum system
    J::Real, # The coupling constant of the spin interaction
    hˣ::Real, # The transverse field
    hᶻ::Real # The longitudinal field
    )::Matrix{ComplexF64}

    𝛔ˣ::Matrix{ComplexF64} = [0 1; 1 0]
    𝛔ᶻ::Matrix{ComplexF64} = [1 0; 0 -1]
    𝒪(args::Pair{Matrix{ComplexF64}, Int}...) = product_local_operators(𝒮, args...)

    d = 𝒮.system_size
    L = 𝒮.local_dimension

    @assert (d == 2) "The local dimension must be 2 in the Ising Model."

    return sum(𝒪(-J/2 * 𝛔ᶻ,i,𝛔ᶻ,i+1) - 𝒪(hˣ*𝛔ˣ + hᶻ*𝛔ᶻ,i) for i in 1:L)
end
export ising_hamiltonian

"""
Generates the matrix corresponding to the local Hamiltonian of the Ising model.

Inputs:
- `𝒮` is the quantum system;
- `ℳ` is an instance of IsingModel;
- `j` is the position of the local Hamiltonian.
"""
function local_hamiltonian(
    𝒮::ExactDiagSystem, # The quantum system
    ℳ::IsingModel, # An instance of IsingModel
    j::Integer, # The position of the local Hamiltonian
    )::Matrix{ComplexF64}

    𝛔ˣ::Matrix{ComplexF64} = [0 1; 1 0]
    𝛔ᶻ::Matrix{ComplexF64} = [1 0; 0 -1]
    𝒪(args::Pair{Matrix{ComplexF64}, Int}...) = product_local_operators(𝒮, args...)

    d = 𝒮.local_dimension
    L = 𝒮.system_size

    @assert (d == 2) "The local dimension must be 2 in the Ising Model."

    j = j ↻ L

    J = ℳ.spin_interaction
    hˣ = ℳ.transverse_field
    hᶻ = ℳ.longitudinal_field

    return 𝒪(-J/4 * 𝛔ᶻ,j-1,𝛔ᶻ,j) - 𝒪(J/4 * 𝛔ᶻ,j,𝛔ᶻ,j+1) - 𝒪(hˣ * 𝛔ˣ + hᶻ * 𝛔ᶻ,j)
end
export local_hamiltonian