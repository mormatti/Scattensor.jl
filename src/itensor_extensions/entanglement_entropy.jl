"""
    entanglement_entropy(mps) -> Float64
    entanglement_entropy(mps, j) -> Vector{Float64}

Computes the von Neumann entanglement entropy of bipartition in a integer position 
j for a ITensor MPS. If `j` is not specified, it computes the entanglement entropy for all sites,
returning a vector of entanglement entropies for each site.

# Example
    julia> sites = siteinds(3,10)
    julia> psi = random_mps(sites)
    julia> entanglement_entropy(psi, 5)
    -1.2813706015259676e-15
"""
function entanglement_entropy(mps::MPS, j::Int)::Float64
    orthogonalize!(mps, j)
    _, S = svd(mps[j], (linkind(mps,j), siteind(mps,j+1)))
    return -sum(p^2 * log(p^2) for p in diag(S)) / log(2)
end

function entanglement_entropy(mps::MPS)::Vector{Float64}
    L = length(mps)
    Sv(j) = entanglement_entropy(mps, j)
    return [Sv(j) for j in 1:L-1]
end

export entanglement_entropy