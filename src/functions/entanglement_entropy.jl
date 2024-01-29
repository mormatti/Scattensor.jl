function entanglement_entropy(ψ::MPS, j::Int)::Float64
    orthogonalize!(ψ, j)
    _,S = svd(ψ[j], (linkind(ψ,j), siteind(ψ,j+1)))
    return -sum(p^2 * log(p^2) for p in diag(S)) / log(2)
end
# export entanglement_entropy later

function entanglement_entropy(ψ::MPS)::Vector{Float64}
    L = length(ψ)
    Sᵥ(j) = entanglement_entropy(ψ, j)
    return [Sᵥ(j) for j in 1:L-1]
end
export entanglement_entropy

