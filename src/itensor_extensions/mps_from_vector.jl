"""
    mps_from_vector(vector, d; kwargs...) -> MPS

Creates a Matrix Product State (MPS) from a vector `v` with local dimension `d`.
The length `L` of the MPS is determined from the size of the vector and `d`.
The keyword arguments `kwargs...` are applied to the `MPS` constructor.

# Example
    julia> v = rand(2^4)
    julia> mps = mps_from_vector(v, 2)
    MPS
    [1] ((dim=2|id=218|"Site,n=1"), (dim=2|id=634|"Link,l=1"))
    [2] ((dim=2|id=634|"Link,l=1"), (dim=2|id=224|"Site,n=2"), (dim=4|id=258|"Link,l=2"))
    [3] ((dim=4|id=258|"Link,l=2"), (dim=2|id=911|"Site,n=3"), (dim=2|id=196|"Link,l=3"))
    [4] ((dim=2|id=196|"Link,l=3"), (dim=2|id=619|"Site,n=4"))
"""
function mps_from_vector(v::Vector, d::Int; cutoff = default_cutoff, maxdim = default_maxdim)::MPS
    # We get L from the vector and d, the check of non-empty vector is done in the function
    L = get_length_from_localdim(v, d)
    # We create the matrix product state (MPS) from the vector v
    sites = siteinds(d, L)
    T = reshape(v, (repeat([d], L)...))
    mps = MPS(ITensor(T, sites...), sites)
    truncate!(mps, cutoff = cutoff, maxdim = maxdim)
    return mps
end

export mps_from_vector