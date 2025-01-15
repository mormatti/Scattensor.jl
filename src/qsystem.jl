mutable struct QuantumSystem
    H0::LocalOperator # The local Hamiltonian
    length::Int # The number of sites of the chain
    localobs::Vector{LocalOperator} # The local observables (Hermitian local operators)
end

function QuantumSystem(H0::LocalOperator, length::Int; localobs::Vector{LocalOperator} = [])
    # It is non sensical if the length of the chain is less than the length of the local Hamiltonian
    if length(H0) > length
        error("The length of the chain must be greater or equal to the length of the local Hamiltonian.")
    end

    return QuantumSystem(H0, length, localobs)
end