"""
    product_locals(L, ops) -> SparseMatrixCSC

Generates a `SparseMatrixCSC` corresponding to the action of a product of local operators, identified by the set `ops`, in a system of length `L`.
The system must be uniform, i.e. all local operators must have the same local dimension `d`.
Each local operator is specified as a pair `(op, j)`, where `op` is the local operator written in the local space (a small matrix) and `j` is the position of the local operator `op`.
Each position is normalized modulo `L`, so that the operator can be applied cyclically.

# Examples
    julia> a = [1 2; 3 4]

    julia> b = [0 1; 1 0]

    julia> product_locals(3,(a,1),(b,3))
    8×8 SparseMatrixCSC{Int64, Int64} with 16 stored entries:
    ⋅  1  ⋅  ⋅  ⋅  2  ⋅  ⋅
    1  ⋅  ⋅  ⋅  2  ⋅  ⋅  ⋅
    ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅  2
    ⋅  ⋅  1  ⋅  ⋅  ⋅  2  ⋅
    ⋅  3  ⋅  ⋅  ⋅  4  ⋅  ⋅
    3  ⋅  ⋅  ⋅  4  ⋅  ⋅  ⋅
    ⋅  ⋅  ⋅  3  ⋅  ⋅  ⋅  4
    ⋅  ⋅  3  ⋅  ⋅  ⋅  4  ⋅
"""
function product_locals(L::Int, ops::Vararg{Tuple{AbstractMatrix, Int}})
    nops = length(ops)
    if nops == 0
        return I  # Identity as a uniform scaling
    end
    # Ensure all matrices match the local space dimension
    d = size(ops[1][1], 1)
    for (op, _) in ops
        if size(op) != (d, d)
            error("In function product_locals each local operator of the list `ops` must have same dimensions. The function is not able to handle different local dimensions.")
        end
    end
    _hilbspace_dimension_warning(SparseMatrixCSC, d, L)  # Ensure the Hilbert space dimension is valid
    n = d^L
    # Normalize positions modulo L and store them
    ops = [(op, pos ↻ L) for (op, pos) in ops]
    # Identity matrix for the local space
    id_matrix = SparseMatrixCSC(Matrix{Int}(I, d, d))
    # Build operators list, filling with identity when no operators are specified
    vectr = [reduce(*, [op for (op, pos) in ops if pos == i]; init = id_matrix) for i in 1:L]
    # Tensor product of the L site operators
    result = reduce(kron, vectr)
    @assert size(result) == (n, n) "The resulting matrix size is incorrect."
    return result
end