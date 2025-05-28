
function operator_reflection(::Type{Matrix}, d::IntTyped, L::IntTypeL) where {IntTypeL <: Integer, IntTyped <: Integer}
    return Matrix(operator_reflection(SparseMatrixCSC, L, d))
end

function operator_reflection(::Type{Hermitian}, d::IntTyped, L::IntTypeL) where {IntTypeL <: Integer, IntTyped <: Integer}
    return Hermitian(operator_reflection(SparseMatrixCSC, L, d))
end

function operator_reflection(::Type{SparseMatrixCSC}, d::IntTyped, L::IntTypeL) where {IntTypeL <: Integer, IntTyped <: Integer}
    # Generate all possible basis states as integer arrays
    basis_states = collect(Iterators.product(ntuple(_ -> 0:d-1, L)...))
    dim_H = d^L  # Total Hilbert space dimension
    R = spzeros(Int, dim_H, dim_H)  # Initialize reflection matrix
    
    # Map basis states to their reflected counterparts
    state_to_index = Dict(state => i for (i, state) in enumerate(basis_states))
    
    for (i, state) in enumerate(basis_states)
        reflected_state = reverse(state)  # Reflect the state
        j = state_to_index[reflected_state]  # Get index of the reflected state
        R[j, i] = 1  # Set matrix element
    end
    
    return R
end

export operator_reflection