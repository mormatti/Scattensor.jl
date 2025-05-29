"""
    apply_reflection(psi::MPS)

Applies the pure reflection of a finite MPS with uniform local dimension.
The pure reflection is defined as the operation that inverts 
the order of the sites in the MPS, effectively reflecting it across the center.

# Example
    julia> sites = siteinds(3,4)
    julia> psi = random_mps(sites)
    julia> phi = apply_reflection(psi)
"""
function apply_reflection(psi::MPS)
    if !is_uniform_localdim(psi)
        error("The input MPS must have have uniform local dimensions.")
    end
    N = length(psi)
    sites = siteinds(psi)
    phi = MPS(N)
    for j in 1:N
        h = N - j + 1
        phi[h] = psi[j] * delta(sites[j], sites[h]')
    end
    noprime!(siteinds, phi)
    truncate!(phi)
    return phi
end

function clean_siteinds!(ψ::MPS)
  new_site_inds = [Index(dim(siteinds(ψ, i)), "Site,n=$i") for i in 1:length(ψ)]
  truncate!(ψ, new_site_inds)
  return ψ
end

export apply_reflection