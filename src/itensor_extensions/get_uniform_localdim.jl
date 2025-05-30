"""
    get_uniform_localdim(siteinds) -> Int
    get_uniform_localdim(mps) -> Int
    get_uniform_localdim(mpo) -> Int

Returns the uniform local dimension of the siteinds of an `MPS` or an `MPO` if all site dimensions are the same.
If the local dimensions are not uniform, it raises an error.

# Examples
    julia> sites = siteinds(3,4)
    julia> get_uniform_localdim(sites)
    3

    julia> dims = [2, 3, 4, 5]
    julia> sites = [Index(d, "Site,n=\$n") for (n, d) in enumerate(dims)]
    julia> psi = random_mps(sites)
    julia> get_uniform_localdim(psi) # Raises an error
"""
function get_uniform_localdim(sites::Vector{<:Index})
    site_dims = [dim(sites[i]) for i in eachindex(sites)]
    all_same = (length(unique(site_dims)) == 1)
    if all_same
        return site_dims[1]  # Return the uniform local dimension
    else
        error("Local dimensions of the MPS not uniform.")
    end
end

function get_uniform_localdim(mps::MPS)
    return get_uniform_localdim(siteinds(mps))
end

function get_uniform_localdim(mpo::MPO)
    d1 = get_uniform_localdim(siteinds_main(mpo))
    d2 = get_uniform_localdim(siteinds_primed(mpo))
    if d1 != d2
        error("Local dimensions of the MPO different between top and bottom indices.")
    end
    return d1  # Return the uniform local dimension
end