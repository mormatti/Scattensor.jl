""" Apply the reflection of a finite MPS with uniform local dimension.
    """
function apply_reflection(ψ::MPS)
    N, Sd = length(ψ), siteinds(ψ)
    ϕ = MPS(N)
    for j in 1:N
        h = N - j + 1
        ϕ[h] = ψ[j] * delta(Sd[j], Sd[h]')
    end
    noprime!(siteinds, ϕ)
    return ϕ
end

export apply_reflection