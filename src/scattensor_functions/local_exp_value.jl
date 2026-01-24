"""
    local_exp_value(A0, v, L0, d, L; pbc=true, addconst=0.0) -> Vector{Float64}

Compute spatially resolved expectation values of a local operator translated along a 1D chain.

Given a local operator `A0` acting on `L0` consecutive sites, this function embeds it into a chain
of length `L` with local dimension `d`, translates it site-by-site using the translation operator,
and computes ⟨v|Aⱼ|v⟩ for each position.

# Arguments
- `A0::Union{Matrix,SparseMatrixCSC}`: Local operator acting on `L0` sites (size `d^L0 × d^L0`).
- `v::Vector{<:Union{Complex,Real}}`: State vector of the full chain (length `d^L`).
- `L0::Integer`: Support size of `A0` in number of sites.
- `d::Integer`: Local Hilbert space dimension.
- `L::Integer`: Total chain length.

# Keyword Arguments
- `pbc`: Periodic boundary condition flag. (Currently not used by the implementation.)
- `addconst::Real=0.0`: Constant added to each returned expectation value (useful for shifting plots).

# Returns
- `Vector{Float64}` of length `L` containing `real(v' * Aⱼ * v) + addconst` for each translated position.

# Notes
- This routine constructs the full embedded operator as a (potentially large) matrix and is intended
  for small systems / exact-vector workflows.
- The `pbc` keyword is currently accepted but not used (see TODO in source).
"""
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
