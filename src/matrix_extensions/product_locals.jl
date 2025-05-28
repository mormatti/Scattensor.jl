# Old code, should be deprecated
""" Generates the whole matrix corresponding to the action of a local operators product.

    ## Inputs
    - `L` is the number of sites of the chain;
    - `d` is the local dimension;
    - `args` is a list of pairs (𝐚, j) where 𝐚 is the local operator written in the local 
    space (small matrix) and j is the position of the local operator 𝐚.
"""
function product_locals(L::Int, d::Int, args::Vararg{Tuple{AbstractMatrix, Int}})
    N = d^L
    
    # Ensure all matrices match the local space dimension
    for (op, _) in args
        @assert size(op) == (d, d) "Each local operator must have dimensions $d x $d."
    end
    
    # Normalize positions modulo L and store them
    args = [(op, pos ↻ L) for (op, pos) in args]
    
    # Identity matrix for the local space
    𝟙 = SparseMatrixCSC(Matrix{Int}(I, d, d))
    
    # Build operators list, filling with identity when no operators are specified
    𝒱 = [reduce(*, [op for (op, pos) in args if pos == i]; init = 𝟙) for i in 1:L]
    
    # Tensor product of the L site operators
    𝐀 = reduce(kron, 𝒱)
    @assert size(𝐀) == (N, N) "The resulting matrix size is incorrect."
    
    return 𝐀
end