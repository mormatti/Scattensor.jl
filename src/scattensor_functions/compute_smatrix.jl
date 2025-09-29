# TODO: all of this code is not yet fully functional, it needs to be tested and debugged.
# TODO: this code still use substitute_siteinds and blind_product_inner... fix it.

function smatrix_real_space(d::Integer, H::MPO, ψ0::MPS, w::MPO, t::Real, N::Integer, χmax::Integer)
    # TODO assert that the local dimension of the MPS is the same as the one of the MPO
    # TODO assert that ψ0 is the groundstate of the Hamiltonian H
    Lψ0 = length(ψ0)
    l = length(w)
    L = Lψ0 - l + 1 # We define the length of the particle positions
    # We construct the single-particle basis wj
    println("Constructing the creation operators for length L...")
    W = Array{MPO}(undef, L)
    for j in 1:L
        W[j] = insert_local(j - 1, deepcopy(w), Lψ0 - (j - 1) - l, d)
        substitute_siteinds!(W[j], ψ0) # FIXME substitution MPS -> MPO siteinds
    end
    println("Constructing the asimptotic out states...")
    ψout = Array{MPS}(undef, L, L)
    for j1 in 1:L
        for j2 in 1:L
            ψtoadd = deepcopy(ψ0)
            ψtoadd = apply(W[j1], ψtoadd; maxdim = χmax)
            ψtoadd = apply(W[j2], ψtoadd; maxdim = χmax)
            normalize!(ψtoadd)
            ψout[j1, j2] = ψtoadd
        end
    end

    println("Constructing the asimptotic in states...")
    ψin = deepcopy(ψout)
    Loffset = Int64(round(Lψ0/4))
    Lmin = 1 + Loffset
    Lmax = L - Loffset


    # We compute the matrix A
    println("Computing the smatrix element expansion in real space...")
    S = zeros(ComplexF64, L, L, L, L) # N+1 because we include the N = 0 term
    for n in 0:N
        Sn = zeros(ComplexF64, L, L, L, L)
        normn = 0
        for j1 in Lmin:Lmax
            for j2 in Lmin:Lmax
                if abs(j1 - j2) < l
                    continue
                end
                factor = 1
                if n > 0
                    factor = -im * t / n
                end
                ψin[j1, j2] = factor * apply(H, ψin[j1, j2]; maxdim = χmax)
                for j1′ in 1:L
                    for j2′ in 1:L
                        if abs(j1′ - j2′) < l
                            continue
                        end
                        term = inner(ψout[j1′, j2′], ψin[j1, j2])
                        S[j1, j2, j1′, j2′] += term
                        Sn[j1, j2, j1′, j2′] = term
                        normn += abs(term)^2
                    end
                end
            end
        end
        println("Norm of the $n-th term: $(sqrt(normn))")
    end
    return S
end

export smatrix_real_space

function smatrix_element_momentum_space(S::Array, k1::Real, k2::Real, k1′::Real, k2′::Real)
    result = 0
    L1 = size(S, 1)
    L2 = size(S, 2)
    L1′ = size(S, 3)
    L2′ = size(S, 4)
    normaliz = sqrt(L1*L2*L1′*L2′)
    for j1 in 1:L1
        for j2 in 1:L2
            for j1′ in 1:L1′
                for j2′ in 1:L2′
                    prefactor = exp(-im * k1 * j1) * exp(-im * k2 * j2)
                    prefactor *= exp(im * k1′* j1′) * exp(im * k2′ * j2′)
                    prefactor /= normaliz
                    result += prefactor * S[j1, j2, j1′, j2′]
                end
            end
        end
    end
    return result
end

export smatrix_momentum_space

function smatrix_real_space_tdvp(
    d::dType, # The local dimension of the system
    H::MPO, # The Hamiltonian of the system
    ψ0::ψ0Type, # The vacuum state
    w::wType, # The Wannier creation operator
    t::tType, # The time at which we want to compute the S-matrix element
    χmax::Int # The maximum bond dimension
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
    L = Lψ0 - l + 1 # We define the length of the particle positions

    # We construct the single-particle basis wj
    println("Constructing the creation operators for length L...")
    W = Array{MPO}(undef, L)
    for j in 1:L
        W[j] = insert_local(j - 1, deepcopy(w), Lψ0 - (j - 1) - l, d)
        substitute_siteinds!(W[j], ψ0)
    end

    println("Constructing the asimptotic out states...")
    ψout = Array{MPS}(undef, L, L)
    for j1 in 1:L
        for j2 in 1:L
            ψtoadd = deepcopy(ψ0)
            ψtoadd = apply(W[j1], ψtoadd; maxdim = χmax)
            ψtoadd = apply(W[j2], ψtoadd; maxdim = χmax)
            normalize!(ψtoadd)
            ψout[j1, j2] = ψtoadd
        end
    end

    println("Constructing the asimptotic in states...")
    ψin = deepcopy(ψout)
    Loffset = Int64(round(Lψ0/4))
    Lmin = 1 + Loffset
    Lmax = L - Loffset

    # We compute the matrix A
    println("Computing the smatrix element expansion in real space...")
    S = zeros(ComplexF64, L, L, L, L) # N+1 because we include the N = 0 term
    for j1 in Lmin:Lmax
        for j2 in Lmin:Lmax
            println("step $j1, $j2")
            if abs(j1 - j2) < l
                continue
            end
            dt = t / 10
            for n in 1:10
                ψin[j1, j2] = tdvp(H, -im * dt, ψin[j1, j2]; maxdim = χmax)
                println("iteration $n, maxlinkdim = $(maxlinkdim(ψin[j1, j2]))")
            end
            for j1′ in 1:L
                for j2′ in 1:L
                    if abs(j1′ - j2′) < l
                        continue
                    end
                    term = inner(ψout[j1′, j2′], ψin[j1, j2])
                    S[j1, j2, j1′, j2′] += term
                end
            end
        end
    end
    return S
end

export smatrix_real_space_tdvp

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

function two_particle_wavefunction(ψ::MPS, ψ0::MPS, W::MPO, d::Int; maxdim::Int = 100)
    normalize!(ψ)
    normalize!(ψ0)
    Lψ = length(ψ)
    L0 = length(ψ0)
    if siteinds(ψ) != siteinds(ψ0)
        substitute_siteinds!(ψ, ψ0)
    end
    @assert Lψ == L0
    Lw = length(W)
    function Wcr(j)
        ret = insert_local(j - 1, deepcopy(W), L0 - (j - 1) - Lw, d)
        substitute_siteinds!(ret, ψ0)
        return ret
    end
    L = Lψ - Lw + 1
    M = zeros(ComplexF64, L, L)
    for j in 1:L
        println("Step: $j")
        ψ0′ = deepcopy(ψ0)
        Wj = Wcr(j)
        ψ0′ = apply(Wj, ψ0′; maxdim = maxdim)
        normalize!(ψ0′)
        for j′ in 1:L
            if abs(j - j′) < Lw
                continue
            end
            ψ0′′ = deepcopy(ψ0′)
            ψ0′′ = apply(Wcr(j′), ψ0′′; maxdim = maxdim)
            normalize!(ψ0′′)
            M[j, j′] = inner(ψ0′′, ψ)
        end
    end
    return M
end
