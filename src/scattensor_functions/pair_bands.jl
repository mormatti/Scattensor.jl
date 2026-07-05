# =============================================================================
# pair_bands.jl — the TWO-PARTICLE sector from bulk matrix elements (B₂)
#
# (MM notes 2026-07-05, both Apple notes.) Basis: pairs of dressed creators at
# ordered positions, |j₁ j₂⟩ = φ̂^{j₁} φ̂^{j₂} |Ω⟩ with j₂ − j₁ ≥ ℓ (disjoint
# supports). The generating data is the four-point function
#
#     H₂(r, r′, m) = ⟨ (c, c+r) | H − E₀ | (c+m, c+m+r′) ⟩     (+ S₂ likewise)
#
# — by bulk translation invariance a function of THREE integers, measurable
# only inside a finite window. The structural fact that makes the spectrum
# computable on an ARBITRARILY LARGE relative domain is cluster decomposition:
#
#     S₂ = s(m)·s(m+r′−r)  +  s(m+r′)·s(m−r)          + connected
#          [direct pairing]   [exchange pairing]
#     H₂ = h(m)s(m+r′−r) + s(m)h(m+r′−r) + (exchange) + connected
#
# where s(r), h(r) are the ONE-body bulk functions (bulk_couplings, B₁) — the
# disconnected skeleton is EXACT in terms of already-measured data, and the
# connected remainder decays exponentially in the distance from the pairing
# manifold (the metric "R" of the note). v1 strategy: tabulate the measured
# window exactly (it contains skeleton + connected), use the skeleton alone
# outside — no basis-function fitting needed; the exponential decay justifies
# the hard crossover.
#
# At fixed total momentum K the center-of-mass Bloch sum gives a generalized
# eigenproblem in the relative coordinate r on a domain of any size Rmax:
# dense two-particle continuum (edges = free kinematics of the measured
# band), bound states (localized eigenvectors), all from one finite chain.
#
# v1: single species (one creator φ). Species generalization: promote the
# s, h vectors and the measured tables to matrices over (α, β) pairs.
# =============================================================================

"""
    pair_couplings(H, Ω, φ, d; seps, hops = 4, center = nothing) -> (S2, H2)

Measure the two-particle bulk data on a finite chain: for a left pair at
`(c, c+r)` and a right pair at `(c+m, c+m+r′)` (vacuum-projected states
`φ̂^{j₁} φ̂^{j₂}|Ω⟩`), return complex arrays indexed `[ri, rj, m + hops + 1]`
with `r = seps[ri]`, `r′ = seps[rj]`, `m ∈ -hops:hops`:

    S2 = ⟨left|right⟩ ,   H2 = ⟨left|(H − E₀)|right⟩ .

`seps` defaults to `ℓ : ℓ+5` (`ℓ` = creator support: minimal non-overlapping
separation). All supports must fit in `[1, L]`; `center` is chosen
automatically mid-chain and can be overridden.
"""
function pair_couplings(H::AbstractMatrix, Ω::AbstractVector, φ::AbstractMatrix, d::Integer;
                        seps = nothing, hops::Integer = 4,
                        center::Union{Nothing, Integer} = nothing)
    L = get_length_from_localdim(size(H, 1), d)
    ℓ = get_length_from_localdim(size(φ, 1), d)
    seps = seps === nothing ? collect(ℓ:(ℓ + 5)) : collect(seps)
    minimum(seps) >= ℓ || error("separations must be ≥ the creator support ℓ = $ℓ (disjoint supports)")
    M = hops
    smax = maximum(seps)
    cmax = L - M - smax - ℓ + 1
    c = center === nothing ? clamp((L - smax - ℓ + 2 - 2M) ÷ 2 + M, M + 1, max(M + 1, cmax)) : Int(center)
    (c - M >= 1 && c + M + smax + ℓ - 1 <= L) ||
        error("pair window does not fit (need L ≥ $(2M + smax + ℓ)): L=$L, center=$c, smax=$smax, hops=$M, support=$ℓ")
    E0 = real(dot(Ω, H * Ω))
    function ket(j1, j2)
        v = embed_operator(φ, d, L, j2) * Ω
        v = embed_operator(φ, d, L, j1) * v
        v -= dot(Ω, v) * Ω
        v
    end
    ns = length(seps)
    lefts = [ket(c, c + r) for r in seps]
    Hlefts = [H * v - E0 * v for v in lefts]
    S2 = zeros(ComplexF64, ns, ns, 2M + 1)
    H2 = zeros(ComplexF64, ns, ns, 2M + 1)
    for (rj, r2) in enumerate(seps), m in -M:M
        v = ket(c + m, c + m + r2)
        for ri in 1:ns
            S2[ri, rj, m + M + 1] = dot(lefts[ri], v)
            H2[ri, rj, m + M + 1] = dot(Hlefts[ri], v)
        end
    end
    S2, H2
end

# one-body function lookup with hard cutoff outside the measured range
_lk(f::AbstractVector, r::Integer, rmax::Integer) = abs(r) <= rmax ? f[r + rmax + 1] : zero(eltype(f))

"""
    pair_blocks(K, s1, h1, S2, H2, seps; hops, Rmax = 300) -> (SK, HK, rgrid)

Assemble the fixed-total-momentum-`K` generalized pair `(S_K, H_K)` on the
relative-coordinate grid `r = minimum(seps) : Rmax`, from:
- the one-body bulk functions `s1, h1` (vectors over `r = -r1b:r1b`, i.e.
  `S[1,1,:]`, `Hc[1,1,:]` of [`bulk_couplings`](@ref)) — used as the exact
  disconnected skeleton (direct + exchange pairings) everywhere;
- the measured tables `S2, H2` of [`pair_couplings`](@ref) — used verbatim
  (skeleton + connected physics) where available.

Center-of-mass convention `|K, r⟩ = Σ_c e^{iK(c + r/2)} |c, c+r⟩`, so
`X_K(r, r′) = Σ_m e^{iK(m + (r′−r)/2)} X₂(r, r′, m)`. Both outputs are
Hermitized. Diagonalize with a Löwdin cutoff (see
[`two_particle_spectrum`](@ref)).
"""
function pair_blocks(K::Real, s1::AbstractVector, h1::AbstractVector,
                     S2::Array{ComplexF64, 3}, H2::Array{ComplexF64, 3},
                     seps; hops::Integer = 4, Rmax::Integer = 300)
    seps = collect(seps)
    smin, smax = extrema(seps)
    r1b = (length(s1) - 1) ÷ 2
    M = hops
    sepidx = Dict(r => i for (i, r) in enumerate(seps))
    rgrid = smin:Rmax
    n = length(rgrid)
    SK = zeros(ComplexF64, n, n)
    HK = zeros(ComplexF64, n, n)
    for (i, r) in enumerate(rgrid), (j, r2) in enumerate(rgrid)
        if haskey(sepidx, r) && haskey(sepidx, r2)
            # measured region: exact window data (skeleton + connected)
            ri, rj = sepidx[r], sepidx[r2]
            for m in -M:M
                ph = exp(im * K * (m + (r2 - r) / 2))
                SK[i, j] += ph * S2[ri, rj, m + M + 1]
                HK[i, j] += ph * H2[ri, rj, m + M + 1]
            end
        else
            # skeleton region: direct + exchange pairings from one-body data
            for m in -(r1b + abs(r2 - r)):(r1b + abs(r2 - r))
                sd1 = _lk(s1, m, r1b);        sd2 = _lk(s1, m + r2 - r, r1b)
                hd1 = _lk(h1, m, r1b);        hd2 = _lk(h1, m + r2 - r, r1b)
                se1 = _lk(s1, m + r2, r1b);   se2 = _lk(s1, m - r, r1b)
                he1 = _lk(h1, m + r2, r1b);   he2 = _lk(h1, m - r, r1b)
                (sd1 == 0 && se1 == 0 && hd1 == 0 && he1 == 0) && continue
                ph = exp(im * K * (m + (r2 - r) / 2))
                SK[i, j] += ph * (sd1 * sd2 + se1 * se2)
                HK[i, j] += ph * (hd1 * sd2 + sd1 * hd2 + he1 * se2 + se1 * he2)
            end
        end
    end
    (SK + SK') / 2, (HK + HK') / 2, rgrid
end

"""
    two_particle_spectrum(K, s1, h1, S2, H2, seps;
                          hops = 4, Rmax = 300, gram_tol = 1e-6,
                          bound_frac = 0.2, bound_weight = 0.7, edge_pad = 1e-3)
        -> (energies, isbound, rgrid)

Diagonalize the fixed-`K` pair problem of [`pair_blocks`](@ref) (Löwdin
cutoff `gram_tol` on the Gram spectrum). Each eigenvector is classified as
BOUND when (i) more than `bound_weight` of its weight (in the orthonormalized
frame) sits within the first `bound_frac` fraction of the relative-distance
grid, AND (ii) its energy lies outside the free-pair kinematic window
`[lo(K) − edge_pad, hi(K) + edge_pad]` computed from the one-body symbols
`ε(k) = h(k)/s(k)`. Criterion (ii) protects against flat-band points (e.g.
XXZ at K = π), where scattering eigenstates are also localized in `r`.
"""
function two_particle_spectrum(K::Real, s1::AbstractVector, h1::AbstractVector,
                               S2::Array{ComplexF64, 3}, H2::Array{ComplexF64, 3},
                               seps; hops::Integer = 4, Rmax::Integer = 300,
                               gram_tol::Real = 1e-6,
                               bound_frac::Real = 0.2, bound_weight::Real = 0.7,
                               edge_pad::Real = 1e-3)
    SK, HK, rgrid = pair_blocks(K, s1, h1, S2, H2, seps; hops, Rmax)
    eS, US = eigen(Hermitian(SK))
    keep = real.(eS) .> gram_tol * maximum(real.(eS))
    W = US[:, keep] * Diagonal(1 ./ sqrt.(real.(eS[keep])))
    Ht = Hermitian(W' * HK * W)
    E, V = eigen(Ht)
    # free-pair kinematic window from the one-body symbols
    r1b = (length(s1) - 1) ÷ 2
    εk(k) = real(sum(exp(-im * k * r) * h1[r + r1b + 1] for r in -r1b:r1b) /
                 sum(exp(-im * k * r) * s1[r + r1b + 1] for r in -r1b:r1b))
    pairE = [εk(K / 2 + q) + εk(K / 2 - q) for q in range(-π, π; length = 601)]
    lo, hi = extrema(pairE)
    ncut = max(1, round(Int, bound_frac * length(rgrid)))
    isbound = falses(length(E))
    for a in eachindex(E)
        (E[a] < lo - edge_pad || E[a] > hi + edge_pad) || continue
        w = W * V[:, a]                        # back to the r-grid frame
        w ./= norm(w)
        isbound[a] = sum(abs2, w[1:ncut]) > bound_weight
    end
    real.(E), isbound, rgrid
end
