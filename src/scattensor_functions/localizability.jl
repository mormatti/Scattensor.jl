# =============================================================================
# localizability.jl — ansatz-free quasiparticle localizability
#
# Reference: MM Kraków notes ("Residues and particle localizability") and the
# Spectroscopy project notes. Core objects, for a normalized state |ψ⟩ and a
# reference |Ω⟩ (usually the ground state) on a chain of L sites with local
# dimension d:
#
#   A(ℓ)        = Tr_{rest} |ψ⟩⟨Ω|          (reduced transition operator,
#                                            support = the last ℓ sites)
#   ‖A‖₁        = max_{‖φ‖∞ ≤ 1} |⟨ψ|φ|Ω⟩|  (orthogonal-Procrustes optimum)
#   C̃(ℓ)       = min(1, ‖A‖₁ √(L/ℓ))       (normalized localizability;
#                                            plateau height = residue,
#                                            C̃(L) = 1 for every state)
#   ℓ_loc       = Σ_{ℓ=0}^{L-1} [1 − C̃(ℓ)]  (mean support of the optimal
#                                            local creator — threshold-free)
#   φ*          = U V† from A = U Σ V†      (the optimal creation operator)
#
# Support convention: the LAST ℓ sites in kron ordering (column-major fastest
# indices). For translation-eigenstates the choice of window is immaterial.
# =============================================================================

"""
    transition_trace_norm(ψ, Ω, ℓ; d=2) -> Float64

Trace norm `‖Tr_rest |ψ⟩⟨Ω|‖₁` of the reduced transition operator on a
support of `ℓ` sites (the last `ℓ` sites in kron ordering).

Efficient at every `ℓ ≤ L−1`: for `ℓ` beyond half chain the outer product is
never formed — with `A = M_ψ M_Ω†` and thin QR `M = Q R`, the singular values
of `A` equal those of `R_ψ R_Ω†` (size `d^(L−ℓ)`).
"""
function transition_trace_norm(ψ::AbstractVector, Ω::AbstractVector, ℓ::Integer; d::Integer = 2)
    L = get_length_from_localdim(length(ψ), d)
    1 <= ℓ <= L - 1 || error("support ℓ must be in 1:$(L-1)")
    m = d^(L - ℓ)
    Mψ = reshape(ψ, d^ℓ, m)
    MΩ = reshape(Ω, d^ℓ, m)
    if d^ℓ <= m
        sum(svdvals(Matrix(Mψ * MΩ')))
    else
        Rψ = qr(Matrix(Mψ)).R
        RΩ = qr(Matrix(MΩ)).R
        sum(svdvals(Rψ * RΩ'))
    end
end

"""
    localizability(ψ, Ω, ℓ; d=2) -> Float64

Normalized localizability `C̃(ℓ) = min(1, ‖Tr_rest|ψ⟩⟨Ω|‖₁ · √(L/ℓ))`.

Rises from 0 and reaches exactly 1 at `ℓ = L` for every normalized state; a
quasiparticle of intrinsic size `ℓ₀` plateaus at its residue for `ℓ ≳ ℓ₀`.
The `√(L/ℓ)` accounts for the coherent harvesting of `ℓ` insertion positions
by the optimal operator (see the Spectroscopy notes).
"""
function localizability(ψ::AbstractVector, Ω::AbstractVector, ℓ::Integer; d::Integer = 2)
    L = get_length_from_localdim(length(ψ), d)
    min(1.0, transition_trace_norm(ψ, Ω, ℓ; d) * sqrt(L / ℓ))
end

"""
    support_size(ψ, Ω; d=2) -> Float64

The threshold-free support size `ℓ_loc = Σ_{ℓ=0}^{L-1} [1 − C̃(ℓ)]` — the
area above the localizability curve, i.e. the mean support of the optimal
local creator of `|ψ⟩` from `|Ω⟩`. Compact quasiparticle ≈ 1–2; weakly bound
composite ≈ its physical size; scattering states ≈ L/2.
"""
function support_size(ψ::AbstractVector, Ω::AbstractVector; d::Integer = 2)
    L = get_length_from_localdim(length(ψ), d)
    s = 1.0                                        # ℓ = 0 contributes 1 − 0
    for ℓ in 1:L-1
        s += 1.0 - localizability(ψ, Ω, ℓ; d)
    end
    s
end

"""
    optimal_creator(ψ, Ω, ℓ; d=2) -> (φstar, residue)

The Procrustes maximizer on the `ℓ`-site support: with
`A = Tr_rest |ψ⟩⟨Ω| = U Σ V†`, the operator `φ* = U V†` (a `d^ℓ × d^ℓ`
unitary) attains `|⟨ψ|φ*|Ω⟩| = ‖A‖₁` — it is the optimal bare creation
operator of `|ψ⟩` on that support (e.g. for quench loading). Also returns
the normalized residue `‖A‖₁ √(L/ℓ)` (unclipped).

Requires `d^ℓ ≤ d^(L−ℓ)` (support at most half the chain), since `φ*` is
materialized as a dense `d^ℓ × d^ℓ` matrix.
"""
function optimal_creator(ψ::AbstractVector, Ω::AbstractVector, ℓ::Integer; d::Integer = 2)
    L = get_length_from_localdim(length(ψ), d)
    m = d^(L - ℓ)
    d^ℓ <= m || error("optimal_creator requires ℓ ≤ L/2 (got ℓ = $ℓ, L = $L)")
    Mψ = reshape(ψ, d^ℓ, m)
    MΩ = reshape(Ω, d^ℓ, m)
    F = svd(Matrix(Mψ * MΩ'))
    F.U * F.Vt, sum(F.S) * sqrt(L / ℓ)
end
