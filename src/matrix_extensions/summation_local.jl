"""
    summation_local(A0, d, L; pbc = false) -> SparseMatrixCSC

Computes the summation of the (p-local) operator `A0` over the chain of length `L`.
Starting from `A0`, constructs the operator `A1 + A2 + A3 + ... + AL` where
`Aj = I ⊗ I ⊗ ... ⊗ A0 ⊗ I ⊗ ... ⊗ I`, I are identities of local dimension d,
and `A0` is placed in the `j`-th position. In other words, it computes the summation ΣⱼAj.
If `pbc` is `true`, the chain is considered periodic, otherwise it is open.

# Example
    julia> sx = [0 1; 1 0]

    julia> summation_local(sx, 2, 4)
    16×16 SparseMatrixCSC{Int64, Int64} with 64 stored entries:
    ⎡⢎⡱⠑⢄⠑⢄⠀⠀⎤
    ⎢⠑⢄⢎⡱⠀⠀⠑⢄⎥
    ⎢⠑⢄⠀⠀⢎⡱⠑⢄⎥
    ⎣⠀⠀⠑⢄⠑⢄⢎⡱⎦
"""
function summation_local(A0::Union{Matrix, SparseMatrixCSC}, d::Integer, L::Integer; pbc::Bool = false)
    # Check if the system size L is a positive integer
    if L <= 0
        error("The system size L must be a positive integer, received: L = $L.")
    end
    # Check if the local dimension d is a positive integer
    if d <= 0
        error("The local dimension d must be a positive integer, received: d = $d.")
    end
    L0 = get_length_from_localdim(A0, d) # We compute the length of the local operator A0
    # Check that L0 is not larger than L
    if L0 > L
        error("Size of the local operator A0 larger than the system size L. Received: L0 = $L0, L = $L.")
    end
    # Define the identity matrix for the complementary space, the translation operator and the Hamiltonian
    Idmext = operator_identity(SparseMatrixCSC, d^(L-L0))
    H0ext = kron(A0, Idmext)
    T = operator_translation(SparseMatrixCSC, d, L) # We obtain the translation operator T
    H = deepcopy(H0ext) # We compute H, the Hamiltonian of the chain of length L
    Lfin = L-1
    if pbc == false
        Lfin = L - L0
    end
    # We loop over the chain and apply the translation operator T to H0ext
    for _ in 1:Lfin
        H0ext = T * H0ext * T'
        H += H0ext
    end
    return H
end