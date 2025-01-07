using Plots
using LinearAlgebra, LinearSolve, SparseArrays
using ITensors, ITensorMPS
using KrylovKit

"""A shortcut binary notation for the periodic modulus."""
function ↻(n::Integer, m::Integer)::Integer
    n > 0 ? (n-1)%m + 1 : m + n%m
end

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

""" Generates the translation operator for a chain of `L` sites with local dimension `d`.
    
    ## Assumptions
     - The system is assumed to be uniform, i.e. the local dimension is the same for all sites.
     - The system, in order to perform a translation, must be in periodic boundary conditions.

    ## Inputs
    - `L` is the number of sites of the chain.
    - `d` is the local dimension.

    ## Outputs
    - The translation operator `T` in matrix form.
    """
function generate_translation_operator_matrix(d::Integer, L::Integer)::Matrix{Int}  
    N::Integer = d^L
    T::Matrix{Int} = zeros(N,N)

    Lst = []
    c = 0
    for _ in 1:d
        lst = []
        for _ in 1:(N/d)
            c = c + 1
            push!(lst, c)
        end
        push!(Lst, lst)
    end

    for indL in eachindex(Lst)
        lst = Lst[indL]
        for ind in eachindex(lst)
            j = lst[ind]
            T[j, ((d*(j-1)+1)%N) + indL - 1] = 1
        end
    end

    return T
end