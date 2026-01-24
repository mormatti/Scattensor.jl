"""
    get_length_from_localdim(N::Integer, d::Integer) -> Integer
    get_length_from_localdim(v::AbstractVector, d::Integer) -> Integer
    get_length_from_localdim(A::AbstractMatrix, d::Integer) -> Integer

Infer the chain length `L` from a total Hilbert space size and a local dimension `d`.

This helper solves `N = d^L` for `L` and validates that the result is (approximately) an integer.
It is useful when you know `d` and you are given a full vector/operator acting on the full Hilbert space.

# Arguments
- `N`: Total Hilbert space dimension.
- `v`: State vector of length `N`.
- `A`: Square operator matrix of size `N × N`.
- `d`: Local dimension per site.

# Returns
- `L::Integer` such that `d^L == N`.

# Examples
    julia> get_length_from_localdim(16, 2)
    4
"""
function get_length_from_localdim(N::Integer, d::Integer)::Integer
    L = log(N) / log(d)
    # Check that L0 is an Integer
    if !(L ≈ round(L))
        error("Size of the local operator matrix incompatible with the local dimension d.")
    end
    # Round L to the nearest integer
    L = Int(round(L))
    return L
end

function get_length_from_localdim(vector::Vector, d::Integer)::Integer
    # Check if vector is a valid vector
    N = length(vector)
    if N == 0
        error("The input vector must not be empty.")
    end
    return get_length_from_localdim(N, d)
end

function get_length_from_localdim(matrix::AbstractMatrix, d::Integer)::Integer
    # Check if matrix is a square matrix
    N = size(matrix)[1]
    N2 = size(matrix)[2]
    if N != N2
        error("The input matrix must be square, i.e. have the same number of rows and columns.")
    end
    if N == 0
        error("The input matrix must not be empty.")
    end
    return get_length_from_localdim(N, d)
end
