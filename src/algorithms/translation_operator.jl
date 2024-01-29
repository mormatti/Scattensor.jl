"""
Generates the translation operator T for a chain of L sites with local dimension d.
The system is assumed to be uniform, i.e. the local dimension is the same for all sites.
The system, in order to perform a translation, must be in periodic boundary conditions.

Inputs:
- `L` is the number of sites of the chain.
- `d` is the local dimension.

Outputs:
- The translation operator `T`.
"""
function translation_operator(
    L::Integer, # The length of the system
    d::Integer # The local dimension
    )::Matrix{ComplexF64}

    N::Integer = d^L
    𝐓::Matrix{ComplexF64} = zeros(N,N)

    Lst = []
    c = 0
    for i ∈ 1:d
        lst = []
        for j in 1:(N/d)
            c = c + 1
            push!(lst, c)
        end
        push!(Lst, lst)
    end

    print(Lst)

    for indL in eachindex(Lst)
        lst = Lst[indL]
        for ind in eachindex(lst)
            j = lst[ind]
            𝐓[j, ((d*(j-1)+1)%N) + indL - 1] = 1
        end
    end

    return 𝐓
end
export translation_operator