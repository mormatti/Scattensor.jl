"""
    operator_identity(AbstractMatrixType, n) -> AbstractMatrix

Returns the identity operator of the specified type `AbstractMatrixType` and size `n x n`, with n `Integer`.
The type can be `Matrix` or `SparseMatrixCSC`.

# Examples
    julia> operator_identity(Matrix, 3)
    3×3 Matrix{Int64}:
    1  0  0
    0  1  0
    0  0  1

    julia> operator_identity(SparseMatrixCSC, 4)
    4×4 SparseMatrixCSC{Int64, Int64} with 4 stored entries:
    1  ⋅  ⋅  ⋅
    ⋅  1  ⋅  ⋅
    ⋅  ⋅  1  ⋅
    ⋅  ⋅  ⋅  1
"""
function operator_identity(::Type{<:Matrix}, n::Integer)::Matrix{Int64}
    _hilbspace_dimension_warning(Matrix, n)
    return Matrix{Int64}(I, n, n)
end

function operator_identity(::Type{<:SparseMatrixCSC}, n::Integer)::SparseMatrixCSC{Int64}
    _hilbspace_dimension_warning(SparseMatrixCSC, n)
    return sparse(I, n, n)
end