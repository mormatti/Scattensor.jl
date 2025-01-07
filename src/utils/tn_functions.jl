"""
    Computes the entanglement entropy in a position j for a ITensor MPS.
    """
function entanglement_entropy(ψ::MPS, j::Int)::Float64
    orthogonalize!(ψ, j)
    _,S = svd(ψ[j], (linkind(ψ,j), siteind(ψ,j+1)))
    return -sum(p^2 * log(p^2) for p in diag(S)) / log(2)
end

"""
    Computes the entanglement entropy array for a ITensor MPS.
    """
function entanglement_entropy(ψ::MPS)::Vector{Float64}
    L = length(ψ)
    Sᵥ(j) = entanglement_entropy(ψ, j)
    return [Sᵥ(j) for j in 1:L-1]
end

"""
    Apply the reflection of a finite MPS with uniform local dimension.
    """
function reflect(ψ::MPS)
    N, Sd = length(ψ), siteinds(ψ)
    ϕ = MPS(N)
    for j in 1:N
        h = N - j + 1
        ϕ[h] = ψ[j] * delta(Sd[j], Sd[h]')
    end
    noprime!(siteinds, ϕ)
    return ϕ
end

"""
    Apply the translation of a finite MPS with uniform local dimension.
    The translation is performed swapping couple of physical indices consecutively.
    """
function translate(ψ::MPS; dir = "right", cutoff = 1e-10)
    L = length(ψ)
    ϕ = copy(ψ)
    if dir == "left"
        for j in 1:L-1
            ϕ = swapbondsites(ϕ, j, cutoff = cutoff)
        end
    elseif dir == "right"
        for j in L-1:-1:1
            ϕ = swapbondsites(ϕ, j, cutoff = cutoff)
        end
    else
        error("Invalid direction")
    end
    return ϕ
end