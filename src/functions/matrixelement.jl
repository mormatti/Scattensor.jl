# Here we define a function to compute the matrix element of a local operator
function matrixelement(ψ1::QuantumState{T}, lo::LocalOperator{T}, j::Inttype, ψ2::QuantumState{T}
    ) where {T <: AbstractMatrix, Inttype <: Int}

    return wavefunction(ψ1) * representation(lo) * wavefunction(ψ2)

end