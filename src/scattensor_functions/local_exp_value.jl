# TODO Check this function and fix the pbc not used (maybe there is a problem there?)
function local_exp_value(
    A0::Union{Matrix, SparseMatrixCSC},
    v::Vector{<:Union{Complex, Real}},
    L0::Integer,
    d::Integer,
    L::Integer;
    pbc = true,
    addconst::Real = 0.0
    )

    @assert size(A0)[1] == size(A0)[2] == d^L0 "Size of A0 incompatible."
    @assert length(v) == d^L "Size of v incompatible."
    @assert L0 <= L "L0 must be less than or equal to L."
    @assert size(v)[1] == d^L "Size of v incompatible."

    A0ext = kron(A0, operator_identity(SparseMatrixCSC, d^(L-L0)))
    T = operator_translation(SparseMatrixCSC, d, L)
    if L%2 == 0
        shift0 = L0/2
    else
        shift0 = (L0-1)/2
    end
    A0ext = (T')^shift0 * A0ext * T^shift0
    vals = [real(v' * A0ext * v) + addconst]
    for _ in 1:L-1
        A0ext = T * A0ext * T'
        push!(vals, real(v' * A0ext * v) + addconst)
    end
    return vals
end

export local_exp_value