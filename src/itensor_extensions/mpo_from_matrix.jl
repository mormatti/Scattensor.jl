"""
    mpo_from_matrix(matrix, d; kwargs...) -> MPO

Creates a `MPO` from a square matrix of size `N x N`, with `N = d^L`, 
where `d` is the local dimension and `L` is the number of sites.
`L` is determined from the size of the matrix and `d`.
The keyword arguments `kwargs...` are applied to the `MPO` constructor.

# Example
    julia> sx = [0 1; 1 0]
    julia> mpo_from_matrix(kron(sx, sx), 2)
    MPO
    [1] ((dim=2|id=746|"Site,n=1"), (dim=2|id=746|"Site,n=1")', (dim=4|id=76|"Link,l=1"))
    [2] ((dim=4|id=76|"Link,l=1"), (dim=2|id=563|"Site,n=2"), (dim=2|id=563|"Site,n=2")')
"""
function mpo_from_matrix(matrix::Matrix, d::Int, cutoff::Real, maxdim::Int)::MPO
    # We get L from the matrix and d.
    # The check of the square matrix size is done in the function
    L = get_length_from_localdim(matrix, d)
    # We create the matrix product operator (MPO) from the matrix A
    sites = siteinds(d, L)
    T = reshape(matrix, (repeat([d], 2L)...))
    mpo = MPO(ITensor(T, (prime.(sites)..., sites...)), sites)
    truncate!(mpo, cutoff = cutoff, maxdim = maxdim)
    return mpo
end

export mpo_from_matrix