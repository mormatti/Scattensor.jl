# TODO write better documentation for this function.
""" Computes the summation of the local operator A0 over the chain of length L.
Starting from A0, constructs the operator A1 + A2 + A3 + ... + AL where
Aj = I ⊗ I ⊗ ... ⊗ A0 ⊗ I ⊗ ... ⊗ I, where I are matrices of local dimension d,
and A0 is in the j-th position.
"""
function summation_local(
    A0::A0Type, 
    L0::L0Type, 
    d::dType, 
    L::LType; 
    pbc::Bool = false
    ) where {
        A0Type <: Union{Hermitian, Matrix, SparseMatrixCSC}, 
        dType <: Integer, 
        L0Type <: Integer,
        LType <: Integer
        }

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