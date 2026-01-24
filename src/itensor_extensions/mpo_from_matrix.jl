"""
    mpo_from_matrix(matrix, d; cutoff=default_cutoff, maxdim=default_maxdim) -> MPO

Create an MPO from a full matrix operator.

The input matrix must be square of size `N × N` with `N = d^L`. The chain length `L` is inferred
from `N` and the local dimension `d`, then the matrix is reshaped into a rank-`2L` tensor and wrapped
as an MPO. Finally, the MPO is truncated.

# Arguments
- `matrix::AbstractMatrix`: Full operator matrix of size `d^L × d^L`.
- `d::Integer`: Local dimension.

# Keyword Arguments
- `cutoff=default_cutoff`: Truncation cutoff used by `truncate!`.
- `maxdim=default_maxdim`: Maximum bond dimension used by `truncate!`.

# Example
    julia> sx = [0 1; 1 0]
    julia> mpo_from_matrix(kron(sx, sx), 2)
    MPO
    [1] ((dim=2|id=746|"Site,n=1"), (dim=2|id=746|"Site,n=1")', (dim=4|id=76|"Link,l=1"))
    [2] ((dim=4|id=76|"Link,l=1"), (dim=2|id=563|"Site,n=2"), (dim=2|id=563|"Site,n=2")')
"""
function mpo_from_matrix(matrix::AbstractMatrix, d::Integer; cutoff::Real = default_cutoff, maxdim::Int = default_maxdim)::MPO
    # We get L from the matrix and d.
    # The check of the square matrix size is done in the function
    L = get_length_from_localdim(matrix, d)
    # We create the matrix product operator (MPO) from the matrix A
    sites = siteinds(d, L)
    T = reshape(Matrix(matrix), (repeat([d], 2L)...))
    mpo = MPO(ITensor(T, (prime.(sites)..., sites...)), sites)
    truncate!(mpo, cutoff = cutoff, maxdim = maxdim)
    return mpo
end
