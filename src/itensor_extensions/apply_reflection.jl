"""
    apply_reflection(mps)

Applies the pure reflection of a finite MPS with uniform local dimension.
The pure reflection is defined as the operation that inverts 
the order of the sites in the MPS, effectively reflecting it across the center.

# Example
    julia> sites = siteinds(3,4)
    julia> psi = random_mps(sites)
    julia> phi = apply_reflection(psi)
"""
function apply_reflection(mps::MPS)
    if !is_uniform_localdim(mps)
        error("The input MPS must have have uniform local dimensions.")
    end
    N = length(mps)
    sites = siteinds(mps)
    phi = MPS(N)
    for j in 1:N
        h = N - j + 1
        phi[h] = mps[j] * delta(sites[j], sites[h]')
    end
    noprime!(siteinds, phi)
    truncate!(phi)
    return phi
end

export apply_reflection