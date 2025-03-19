# Old code, should be deprecated
""" Generates the whole matrix corresponding to the action of a local operators product.

    ## Inputs
    - `L` is the number of sites of the chain;
    - `d` is the local dimension;
    - `args` is a list of pairs (𝐚, j) where 𝐚 is the local operator written in the local 
    space (small matrix) and j is the position of the local operator 𝐚.
    """    
function _product_locals(L::Int, d::Int, args::Vararg{Tuple{AbstractMatrix, Int}})
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

function product_direct(M1::TrainType, M2::TrainType) where {TrainType <: Union{MPS, MPO}}
    L1, L2 = length(M1), length(M2) # Lengths of M1 and M2
    M = MPO(L1 + L2) # Combined MPO
    for i in 1:L1 # Copy tensors from M1
        M[i] = M1[i]
    end
    for j in 1:L2 # Copy tensors from M2
        M[L1 + j] = M2[j]
    end
    # Set the link between M[L1] and M[L1 + 1] to a trivial scalar index (dimension 1)
    trivial_link = Index(1, "Link,trivial")
    v1 = ITensor(trivial_link)
    v2 = ITensor(trivial_link)
    v1[trivial_link=>1] = 1
    v2[trivial_link=>1] = 1
    M[L1] = M[L1] * v1 # Update M[L1] to connect to trivial link
    M[L1 + 1] = M[L1 + 1] * v2 # Update M[L1 + 1] to connect to trivial link
    return M
end

export product_direct

function product_inner(mps1::MPS, mps2::MPS)
    return inner(mps1', replace_siteinds(mps2, siteinds(mps1)))
end

export product_inner

function product_outer(mps1::MPS, mps2::MPS)
    return outer(mps1', replace_siteinds(mps2, siteinds(mps1)))
end

export product_outer

function product_matricial(mpo::MPO, mps::MPS)
    substitute_siteinds!(mpo, mps)
    return apply(mpo, mps)
end

function product_matrcial(mpo1::MPO, mpo2::MPS)
    substitute_siteinds!(mpo1, mpo2)
    return apply(mpo1, mpo2)
end

export product_matricial

function kron_power(A, n)
    result = A
    for _ in 2:n
        result = kron(result, A)
    end
    return result
end

export kron_power