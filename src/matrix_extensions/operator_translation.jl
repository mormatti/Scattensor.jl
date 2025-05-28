# Old version, deprecated
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
function _generate_translation_operator_matrix(d::Integer, L::Integer)::SparseMatrixCSC{Int64, Int64}  
    N = d^L
    T = spzeros(Float64, N, N)

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

function operator_translation(::Type{Matrix}, d::IntTyped, L::IntTypeL) where {IntTypeL <: Integer, IntTyped <: Integer}
    return Matrix(operator_translation(SparseMatrixCSC, L, d))
end

function operator_translation(::Type{SparseMatrixCSC}, d::IntTyped, L::IntTypeL) where {IntTypeL <: Integer, IntTyped <: Integer}
    # Generate all possible basis states as integer arrays
    basis_states = collect(Iterators.product(ntuple(_ -> 0:d-1, L)...))
    dim_H = d^L  # Total Hilbert space dimension
    T = spzeros(Int, dim_H, dim_H)  # Initialize reflection matrix
    
    # Map basis states to their reflected counterparts
    state_to_index = Dict(state => i for (i, state) in enumerate(basis_states))
    
    for (i, state) in enumerate(basis_states)
        reflected_state = Tuple(circshift([i for i in state], 1))  # Reflect the state
        j = state_to_index[reflected_state]  # Get index of the reflected state
        T[i, j] = 1  # Note the definition: T|j⟩ = |j+1⟩
    end
    
    return T
end