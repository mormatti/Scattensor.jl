
"""
    operator_reflection(AbstractMatrixType, d, L) -> AbstractMatrix

Returns the reflection operator of the specified type `AbstractMatrixType` for a Hilbert space of dimension `d^L`, where `d` is an `Integer` representing the local dimension of each subsystem and `L` is an `Integer` representing the number of subsystems.
The reflection operator is computed in the canonical basis, i.e. in the tensor product basis of the subsystems in the choosen order.

# Examples
    julia> operator_reflection(Matrix, 2, 2)
    4×4 Matrix{Int64}:
    1  0  0  0
    0  0  1  0
    0  1  0  0
    0  0  0  1

    julia> operator_reflection(SparseMatrixCSC, 2, 3)
    8×8 SparseMatrixCSC{Int64, Int64} with 8 stored entries:
    1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
    ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅
    ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅
    ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅
    ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
    ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅
    ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅
    ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1
"""
function operator_reflection(::Type{<:Matrix}, d::Integer, L::Integer)::Matrix{Int64}
    _hilbspace_dimension_warning(Matrix, d, L)
    return Matrix(operator_reflection(SparseMatrixCSC, d, L)) # delegate to SparseMatrixCSC implementation
end

function operator_reflection(::Type{<:SparseMatrixCSC}, d::Integer, L::Integer)::SparseMatrixCSC{Int64}
    _hilbspace_dimension_warning(SparseMatrixCSC, d, L)
    # Generate all possible basis states as integer arrays
    basis_states = collect(Iterators.product(ntuple(_ -> 0:d-1, L)...))
    dim_H = d^L  # Total Hilbert space dimension
    R = spzeros(Int, dim_H, dim_H)  # Initialize reflection matrix
    # Map basis states to their reflected counterparts
    state_to_index = Dict(state => i for (i, state) in enumerate(basis_states))
    # Iterate over each basis state and fill the reflection matrix
    for (i, state) in enumerate(basis_states)
        reflected_state = reverse(state)  # Reflect the state
        j = state_to_index[reflected_state]  # Get index of the reflected state
        R[j, i] = 1  # Set matrix element
    end
    return R
end

export operator_reflection