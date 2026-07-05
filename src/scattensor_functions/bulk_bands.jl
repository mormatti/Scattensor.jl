# =============================================================================
# bulk_bands.jl — full dispersion relations from BULK matrix elements
#
# (MM note 2026-07-05.) Given dressed creation operators φ̂_α extracted once
# (e.g. from trap seeds via `optimal_creator`), translate the OPERATORS along
# the bulk of a (possibly open) chain: B₁ = { φ̂_α^j |Ω⟩ }. Translation
# invariance of the bulk makes overlaps and Hamiltonian elements functions of
# the separation only,
#
#     S_αβ(r) = ⟨φ_α^i|φ_β^{i+r}⟩ ,   H_αβ(r) = ⟨φ_α^i|(H−E0)|φ_β^{i+r}⟩ ,
#
# both decaying exponentially in |r| for gapped systems. The Fourier symbols
# S(k), H(k) then give the CONTINUOUS bands from the generalized eigenproblem
#
#     H(k) c = ε(k) S(k) c        (Wannier interpolation / CORE-style),
#
# on any momentum grid — no 2π/L sampling. Works in OBC (no T needed).
# =============================================================================

"""
    embed_operator(op, d, L, j) -> SparseMatrixCSC

Place a local operator `op` (acting on `ℓ = log_d(size(op,1))` contiguous
sites) at position `j` (1-based, leftmost site of the support) of a chain of
`L` sites: `1_{j-1} ⊗ op ⊗ 1_{L-j-ℓ+1}` in kron ordering.
"""
function embed_operator(op::AbstractMatrix, d::Integer, L::Integer, j::Integer)
    ℓ = get_length_from_localdim(size(op, 1), d)
    1 <= j <= L - ℓ + 1 || error("operator of support $ℓ does not fit at site $j of $L")
    left = operator_identity(SparseMatrixCSC, d^(j - 1))
    right = operator_identity(SparseMatrixCSC, d^(L - j - ℓ + 1))
    kron(left, kron(sparse(op), right))
end

"""
    bulk_couplings(H, Ω, creators, d; center, rmax) -> (S, Hc)

Bulk matrix elements between translated dressed states
`|α, j⟩ = φ̂_α^j |Ω⟩` (vacuum-connected: the component along `|Ω⟩` is
projected out). Returns two `N × N × (2 rmax + 1)` arrays indexed by
`[α, β, r + rmax + 1]` for separations `r = -rmax … rmax`:

    S[α,β,·]  = ⟨α, c | β, c+r⟩
    Hc[α,β,·] = ⟨α, c | (H − E0) | β, c+r⟩

with the left state pinned at `center` (defaults to mid-chain). Choose
`center` and `rmax` so all supports stay in the translation-invariant bulk;
elements decay exponentially in `|r|` for gapped `H`.
"""
function bulk_couplings(H::AbstractMatrix, Ω::AbstractVector, creators::Vector, d::Integer;
                        center::Union{Nothing, Integer} = nothing, rmax::Integer = 4)
    L = get_length_from_localdim(size(H, 1), d)
    ℓs = [get_length_from_localdim(size(op, 1), d) for op in creators]
    c = center === nothing ? (L - maximum(ℓs)) ÷ 2 + 1 : Int(center)
    (c - rmax >= 1 && c + rmax + maximum(ℓs) - 1 <= L) ||
        error("window [center ± rmax + support] exceeds the chain; reduce rmax or enlarge L")
    E0 = real(dot(Ω, H * Ω))
    N = length(creators)
    # dressed states at every needed position, vacuum component removed
    ket(α, j) = begin
        v = embed_operator(creators[α], d, L, j) * Ω
        v -= dot(Ω, v) * Ω
        v
    end
    S = zeros(ComplexF64, N, N, 2rmax + 1)
    Hc = zeros(ComplexF64, N, N, 2rmax + 1)
    lefts = [ket(α, c) for α in 1:N]
    Hlefts = [H * lv - E0 * lv for lv in lefts]
    for r in -rmax:rmax, β in 1:N
        v = ket(β, c + r)
        for α in 1:N
            S[α, β, r + rmax + 1] = dot(lefts[α], v)
            Hc[α, β, r + rmax + 1] = dot(Hlefts[α], v)
        end
    end
    S, Hc
end

"""
    wannier_bands(S, Hc; nk = 200, gram_tol = 1e-8) -> (ks, bands)

Continuous band structure from bulk couplings (output of
[`bulk_couplings`](@ref)): build the Fourier symbols
`S(k) = Σ_r e^{-ikr} S(r)`, `H(k) = Σ_r e^{-ikr} Hc(r)` and solve the
generalized eigenproblem `H(k) c = ε S(k) c` on `nk` momenta in `[-π, π)`.
Near-singular directions of `S(k)` (Gram eigenvalues below
`gram_tol · max`) are discarded (Löwdin cutoff). Returns the momentum grid
and an `nk × N` matrix of bands (NaN where a band is cut).
"""
function wannier_bands(S::Array{ComplexF64, 3}, Hc::Array{ComplexF64, 3};
                       nk::Integer = 200, gram_tol::Real = 1e-8)
    N = size(S, 1)
    rmax = (size(S, 3) - 1) ÷ 2
    ks = collect(range(-π, π; length = nk + 1))[1:end-1]
    bands = fill(NaN, nk, N)
    for (ik, k) in enumerate(ks)
        Sk = zeros(ComplexF64, N, N); Hk = zeros(ComplexF64, N, N)
        for r in -rmax:rmax
            ph = exp(-im * k * r)
            @views Sk .+= ph .* S[:, :, r + rmax + 1]
            @views Hk .+= ph .* Hc[:, :, r + rmax + 1]
        end
        Sk = (Sk + Sk') / 2; Hk = (Hk + Hk') / 2
        eS, US = eigen(Hermitian(Sk))
        keep = real.(eS) .> gram_tol * maximum(real.(eS))
        any(keep) || continue
        W = US[:, keep] * Diagonal(1 ./ sqrt.(real.(eS[keep])))
        ε = eigvals(Hermitian(W' * Hk * W))
        bands[ik, 1:length(ε)] .= real.(ε)
    end
    ks, bands
end

"""
    energy_variance(H, ψ) -> Float64

`⟨H²⟩ − ⟨H⟩²` of a normalized state: the certificate that an effective-basis
eigenstate is a true eigenstate (variance → 0 exponentially with the basis
range) rather than a continuum admixture (variance stays O(1)).
"""
function energy_variance(H::AbstractMatrix, ψ::AbstractVector)
    Hψ = H * ψ
    real(dot(Hψ, Hψ)) - real(dot(ψ, Hψ))^2
end
