"""
    operator_translation(AbstractMatrixType, d, L) -> AbstractMatrix

Returns the translation operator of the specified type `AbstractMatrixType` for a Hilbert space of dimension `d^L`, where `d` is an `Integer` representing the local dimension of each subsystem and `L` is an `Integer` representing the number of subsystems.
The translation operator is computed in the canonical basis, i.e. in the tensor product basis of the subsystems in the choosen order.

# Examples
    julia> operator_translation(Matrix, 2, 2)
    4×4 Matrix{Int64}:
    1  0  0  0
    0  0  1  0
    0  1  0  0
    0  0  0  1

    julia> operator_translation(SparseMatrixCSC, 3, 2)
    9×9 SparseMatrixCSC{Int64, Int64} with 9 stored entries:
    1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
    ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅
    ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅
    ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
    ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅
    ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅
    ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
    ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅
    ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1
"""
function operator_translation(::Type{Matrix}, d::Integer, L::Integer)
    _hilbspace_dimension_warning(Matrix, d, L)
    return Matrix(operator_translation(SparseMatrixCSC, d, L)) # delegate to SparseMatrixCSC implementation
end

function operator_translation(::Type{SparseMatrixCSC}, d::Integer, L::Integer)
    _hilbspace_dimension_warning(SparseMatrixCSC, d, L)
    # Generate all possible basis states as integer arrays
    basis_states = collect(Iterators.product(ntuple(_ -> 0:d-1, L)...))
    dim_H = d^L  # Total Hilbert space dimension
    T = spzeros(Int, dim_H, dim_H)  # Initialize translation matrix
    # Map basis states to their reflected counterparts
    state_to_index = Dict(state => i for (i, state) in enumerate(basis_states))
    # Iterate over each basis state and fill the translation matrix
    for (i, state) in enumerate(basis_states)
        reflected_state = Tuple(circshift([i for i in state], 1))  # Reflect the state
        j = state_to_index[reflected_state]  # Get index of the reflected state
        T[i, j] = 1  # Note the definition: T|j⟩ = |j+1⟩
    end
    return T
end