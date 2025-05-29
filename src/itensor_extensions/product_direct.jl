"""
    kron(mps1::MPS, mps2::MPS) -> MPS
    kron(mpo1::MPO, mpo2::MPO) -> MPO

    Kronecker product of two MPS or MPO objects.
    The resulting MPS or MPO has the length equal to the sum of the lengths of the two inputs.
    The link between the last tensor of the first MPS/MPO and the first tensor of the second MPS/MPO is set to a trivial scalar index (dimension 1).
"""
function Base.kron(M1::TrainType, M2::TrainType) where {TrainType <: Union{MPS, MPO}}
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

export kron

kron