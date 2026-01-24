"""
    entanglement_entropy(mps) -> Float64
    entanglement_entropy(mps, j) -> Float64

Compute the von Neumann entanglement entropy across a bipartition of an MPS.

If `j` is provided, computes the entropy of the cut between sites `j` and `j+1`.
If `j` is not provided, computes the entropy for all cuts and returns a vector of length `L-1`.

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
