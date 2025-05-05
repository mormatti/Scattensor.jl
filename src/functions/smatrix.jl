function smatrix_expansion_real_space(
    d::dType, # The local dimension of the system
    H::MPO, # The Hamiltonian of the system
    ψ0::ψ0Type, # The vacuum state
    w::wType, # The Wannier creation operator
    N::NType, # The number of terms of the exponential expansion
    χmax::Int # The maximum bond dimension
    ) where {
        dType <: Integer,
        ψ0Type <: MPS,
        wType <: MPO,
        NType <: Int
        }

    # TODO assert that the local dimension of the MPS is the same as the one of the MPO
    # TODO assert that ψ0 is the groundstate of the Hamiltonian H

    Lψ0 = length(ψ0)
    l = length(w)
    L = Int64(Lψ0/2) - l + 1 # We define the length of the particle positions

    # We construct the single-particle basis wj
    println("Constructing the single-particle basis wj...")
    W = Array{MPO}(undef, L)
    Loffset = Int64(Lψ0/4)
    for j in 1:L
        W[j] = insert_local(j - 1 + Loffset, deepcopy(w), Lψ0 - Loffset - (j - 1) - l, d)
        substitute_siteinds!(W[j], ψ0)
    end

    # We construct the asimptotic in states
    println("Constructing the asimptotic in states...")
    ψin = Array{MPS}(undef, L, L)
    for j1 in 1:L
        for j2 in 1:L
            ψtoadd = deepcopy(ψ0)
            ψtoadd = apply(W[j1], ψtoadd; maxdim = χmax)
            ψtoadd = apply(W[j2], ψtoadd; maxdim = χmax)
            normalize!(ψtoadd)
            ψin[j1, j2] = ψtoadd
        end
    end

    # We construct the asimptotic out states as a deepcopy of the in states
    # We need also to substitute the site indices of the out states
    # println("Constructing the asimptotic out states...")
    ψout = deepcopy(ψin)
    # It seems that the following code is not needed...
    # for j1 in 1:L
    #     for j2 in 1:L
    #         substitute_siteinds!(ψout[j1, j2], ψ0)
    #     end
    # end

    # We compute the matrix A
    println("Computing the smatrix element expansion in real space...")
    A = zeros(ComplexF64, L, L, L, L, N+1) # N+1 because we include the N = 0 term
    println("Total steps = $N, Current step: ")
    for n in 0:N
        print(",", n)
        for j1 in 1:L
            for j2 in 1:L
                # We ensure that the two j indices are not too close
                if abs(j1 - j2) < l
                    continue
                end
                ψin[j1, j2] = apply(H, ψin[j1, j2]; maxdim = χmax)
                for j1′ in 1:L
                    for j2′ in 1:L
                        if abs(j1′ - j2′) < l
                            continue
                        end
                        A[j1, j2, j1′, j2′, n+1] = inner(ψout[j1′, j2′], ψin[j1, j2])
                    end
                end
            end
        end
    end

    return A
end

export smatrix_expansion_real_space

# Given A with 4 j indices and 1 n index, we can compute an S-matrix element in the momentum space
function smatrix_element_momentum_space(A, k1, k2, k1′, k2′, t; N::Int = -1)
    result = 0
    L = size(A, 1)
    if N < 0
        N = size(A, 5) - 1 # We need to subtract 1 because we included the N = 0 term
    end
    for n in 0:N
        for j1 in 1:L
            for j2 in 1:L
                for j1′ in 1:L
                    for j2′ in 1:L
                        result += (-im * t)^n / (factorial(n)) * exp(-im * k1 * j1) * exp(-im * k2 * j2) * exp(im * k1′* j1′) * exp(im * k2′ * j2′) * A[j1, j2, j1′, j2′, n+1]
                    end
                end
            end
        end
    end
    return result
end

export smatrix_element_momentum_space

function smatrix_expansion_real_space_tdvp(
    d::dType, # The local dimension of the system
    H::MPO, # The Hamiltonian of the system
    ψ0::ψ0Type, # The vacuum state
    w::wType, # The Wannier creation operator
    χmax::Int, # The maximum bond dimension
    t::tType # The time at which we want to compute the S-matrix element
    ) where {
        dType <: Integer,
        ψ0Type <: MPS,
        wType <: MPO,
        tType <: Real
        }

    # TODO assert that the local dimension of the MPS is the same as the one of the MPO
    # TODO assert that ψ0 is the groundstate of the Hamiltonian H

    Lψ0 = length(ψ0)
    l = length(w)
    L = Int64(Lψ0/2) - l + 1 # We define the length of the particle positions

    # We construct the single-particle basis wj
    println("Constructing the single-particle basis wj...")
    W = Array{MPO}(undef, L)
    Loffset = Int64(Lψ0/4)
    for j in 1:L
        W[j] = insert_local(j - 1 + Loffset, deepcopy(w), Lψ0 - Loffset - (j - 1) - l, d)
        substitute_siteinds!(W[j], ψ0)
    end

    # We construct the asimptotic in states
    println("Constructing the asimptotic in states...")
    ψin = Array{MPS}(undef, L, L)
    for j1 in 1:L
        for j2 in 1:L
            ψtoadd = deepcopy(ψ0)
            ψtoadd = apply(W[j1], ψtoadd; maxdim = χmax)
            ψtoadd = apply(W[j2], ψtoadd; maxdim = χmax)
            normalize!(ψtoadd)
            ψin[j1, j2] = ψtoadd
        end
    end

    # We construct the asimptotic out states as a deepcopy of the in states
    # We need also to substitute the site indices of the out states
    # println("Constructing the asimptotic out states...")
    ψout = deepcopy(ψin)
    # It seems that the following code is not needed...
    # for j1 in 1:L
    #     for j2 in 1:L
    #         substitute_siteinds!(ψout[j1, j2], ψ0)
    #     end
    # end

    # We compute the matrix A
    println("Computing the smatrix element expansion in real space...")
    A = zeros(ComplexF64, L, L, L, L) # N+1 because we include the N = 0 term
    for j1 in 1:L
        for j2 in 1:L
            if abs(j1 - j2) < l
                continue
            end
            println("TDVP Step: $j1, $j2")
            ψin[j1, j2] = tdvp(H, -im * t, ψin[j1, j2]; maxdim = χmax)
            for j1′ in 1:L
                for j2′ in 1:L
                    if abs(j1′ - j2′) < l
                        continue
                    end
                    A[j1, j2, j1′, j2′] = inner(ψout[j1′, j2′], ψin[j1, j2])
                end
            end
        end
    end

    return A
end

export smatrix_expansion_real_space_tdvp

# Given A with 4 j indices and 1 n index, we can compute an S-matrix element in the momentum space
function smatrix_element_momentum_space_tdvp(A, k1, k2, k1′, k2′)
    result = 0
    L = size(A, 1)
    for j1 in 1:L
        for j2 in 1:L
            for j1′ in 1:L
                for j2′ in 1:L
                    result += exp(-im * k1 * j1) * exp(-im * k2 * j2) * exp(im * k1′* j1′) * exp(im * k2′ * j2′) * A[j1, j2, j1′, j2′]
                end
            end
        end
    end
    return result
end

export smatrix_element_momentum_space_tdvp