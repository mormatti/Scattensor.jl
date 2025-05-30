"""
    is_uniform_localdim(siteinds) -> Bool
    is_uniform_localdim(mps) -> Bool
    is_uniform_localdim(mpo) -> Bool

Returns true if the uniform local dimension of the siteinds of an `MPS` or an `MPO` if all site dimensions are the same.

# Examples
    julia> sites = siteinds(3,4)
    julia> is_uniform_localdim(sites)
    true

    julia> dims = [2, 3, 4, 5]
    julia> sites = [Index(d, "Site,n=\$n") for (n, d) in enumerate(dims)]
    julia> psi = random_mps(sites)
    julia> is_uniform_localdim(psi)
    false
"""
function is_uniform_localdim(sites::Vector{<:Index})::Bool
    site_dims = [dim(sites[i]) for i in eachindex(sites)]
    all_same = (length(unique(site_dims)) == 1)
    if all_same
        return true
    else
        return false
    end
end

function is_uniform_localdim(mps::MPS)::Bool
    return is_uniform_localdim(siteinds(mps))
end

function is_uniform_localdim(mpo::MPO)::Bool
    return is_uniform_localdim(siteinds(mpo))
end

function is_uniform_localdim(mpo::MPO)
    d1 = get_uniform_localdim(siteinds_main(mpo))
    d2 = get_uniform_localdim(siteinds_primed(mpo))
    if d1 != d2
        return false
    end
    return true
end