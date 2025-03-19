function substitute_siteinds!(mpo1::MPO, mpo2::MPO)
    @assert length(mpo1) == length(mpo2) "The MPOs must have the same length."
    si1 = siteinds(mpo1)
    si2 = siteinds(mpo2)
    for j in eachindex(mpo1)
        mpo1[j] = mpo1[j] * delta(si1[j][1], si2[j][1])
        mpo1[j] = mpo1[j] * delta(si1[j][2], si2[j][2])
        # TODO: nonsense to do that (it affects performance), one should do sometinhg like:
        # ITensors.replaceind(mpo1[j], si1[j][1], si2[j][1])
        # ITensors.replaceind(mpo1[j], si1[j][2], si2[j][2])
        # But dunno it does not work! Fix it.
    end
end

function substitute_siteinds!(mps1::MPS, mps2::MPS)
    @assert length(mps1) == length(mps2) "The two MPS must have the same length."
    si1 = siteinds(mps1)
    si2 = siteinds(mps2)
    for j in eachindex(mps1)
        mps1[j] = mps1[j] * delta(si1[j], si2[j])
    end
end

function substitute_siteinds!(mps::MPS, mpo::MPO)
    error("Don't use this method, substitute mpo indices with the mps one instead.")
end

function substitute_siteinds!(mpo::MPO, mps::MPS)
    @assert length(mpo) == length(mps) "The MPS and the MPO must have the same length."
    simpo = siteinds(mpo)
    simps = siteinds(mps)
    for j in eachindex(mpo)
        mpo[j] = mpo[j] * delta(simpo[j][1], prime(simps[j]))
        mpo[j] = mpo[j] * delta(simpo[j][2], simps[j])
    end
end

export substitute_siteinds!

function insert_local(L1::L1Type, localmpo::MPO, L2::L2Type, d::dType) where {L1Type <: Integer, L2Type <: Integer, dType <: Integer}
    idm(L) = operator_identity(MPO, d, L)
    if L1 == 0 && L2 == 0
        return c
    elseif L1 == 0
        return product_direct(localmpo, idm(L2))
    elseif L2 == 0
        return product_direct(idm(L1), localmpo)
    else
        return product_direct(idm(L1), product_direct(localmpo, idm(L2)))
    end
end 

export insert_local

# TODO write better documentation for this function.
""" Computes the summation of the local operator A0 over the chain of length L.
Starting from A0, constructs the operator A1 + A2 + A3 + ... + AL where
Aj = I ⊗ I ⊗ ... ⊗ A0 ⊗ I ⊗ ... ⊗ I, where I are matrices of local dimension d,
and A0 is in the j-th position.
"""
function summation_local(A0::A0Type, L0::L0Type, d::dType, L::LType; pbc::Bool = false) where {A0Type <: Union{Hermitian, Matrix, SparseMatrixCSC}, dType <: Integer, L0Type <: Integer, LType <: Integer}
    @assert size(A0)[1] == size(A0)[2] == d^L0 "Size of A0 incompatible."
    Idmext = operator_identity(SparseMatrixCSC, d^(L-L0)) # Def: the identity matrix for the complementary space
    H0ext = kron(A0, Idmext)
    T = operator_translation(SparseMatrixCSC, d, L) # We obtain the translation operator T
    H = deepcopy(H0ext) # We compute H, the Hamiltonian of the chain of length L
    Lfin = L-1
    if pbc == false
        Lfin = L - L0
    end
    for _ in 1:Lfin
        H0ext = T * H0ext * T'
        H += H0ext
    end
    if A0Type <: Hermitian
        H = Hermitian(H)
    end
    return H
end

function summation_local(mpo::MPO, d::dType, L::LType; func::Function = identity, cutoff = 10^(-12)) where {LType <: Integer, dType <: Integer}
    L0 = length(mpo)
    Lc = L - L0
    if func == identity
        f(j) = 1
        func = f
    end
    finalsum = func(0) * insert_local(0, mpo, Lc, d)
    for j in 1:Lc
        tosum = func(j) * insert_local(j, mpo, Lc - j, d)
        substitute_siteinds!(tosum, finalsum)
        finalsum += tosum
    end
    finalsum = truncate(finalsum, cutoff = cutoff)
end

export summation_local

function local_expvals(ψ::MPS, A₀::MPO, d::dType) where {dType <: Integer}
    if length(A₀) < length(ψ)
        L₀ = length(A₀)
        L = length(ψ)
        Lc = L - L₀
        vals = []
        for j in 1:Lc
            Aext = insert_local(j, A₀, Lc - j, d)
            substitute_siteinds!(Aext, ψ)
            push!(vals, real(inner(ψ', Aext, ψ)))
        end
        return vals
    elseif length(A₀) == length(ψ)
        return real(inner(mps', A₀, mps))
    else
        error("The length of the local operator cannot be greater than the one of the .")
    end
end

function local_expvals(ψ::Vector{MPS}, A₀::MPO, d::dType) where {dType <: Integer}
    # We assert that all the MPS inside ψ has the same length
    @assert all(mps -> length(mps) == length(first(ψ)), ψ) "All the MPS must have the same length"
    
    matrix = zeros(length(ψ), length(ψ[1]) - length(A₀))
    for α in eachindex(ψ)
        matrix[α,:] = local_expvals(ψ[α], A₀, d)
    end
    return matrix
end

export local_expvals