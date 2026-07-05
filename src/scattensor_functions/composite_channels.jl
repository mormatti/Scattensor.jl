# =============================================================================
# composite_channels.jl — COMPOSITE CHANNELS in the two-particle sector and the
# ab-initio FRIEDRICHS–LEE reduction
#
# A composite channel is an extra one-body-LIKE channel added to the fixed-K
# pair problem of pair_bands.jl: a localized multi-site creator χ̂ (e.g. a
# 5-site window operator that creates a bound pair directly, with overlapping
# constituents) whose Bloch sum
#
#     |K, χ⟩ = Σ_c e^{iKc} χ̂_c |Ω⟩
#
# is appended to the relative-coordinate pair basis
#
#     |K, r⟩ = Σ_c e^{iK(c + r/2)} φ̂_c φ̂_{c+r} |Ω⟩ ,   r ≥ ℓ (disjoint supports).
#
# The generalized eigenproblem (the overlap S includes all cross terms) handles
# the non-orthogonality — and even exact linear dependence — of the enlarged
# basis: Löwdin null directions are simply dropped. Two payoffs:
#
#  (1) VARIATIONAL: the composite fills the short-distance region r < ℓ that
#      the disjoint-support pair basis cannot reach (overlapping constituents,
#      genuine bound-pair structure).
#
#  (2) STRUCTURAL: orthogonalizing the composite against the pair continuum
#      yields the Friedrichs–Lee data ab initio — a bare channel energy
#      E_χ(K) = ⟨χ̃|H|χ̃⟩ and couplings g_a = ⟨χ̃|H|a⟩ to the discrete
#      continuum levels ε_a — from which the resonance position, its width
#      Γ(K), and a SMOOTH resonant phase shift δ(E) follow without any
#      finite-volume level-crossing analysis.
#
# Conventions follow pair_bands.jl exactly: center-of-mass phase
# e^{iK(m + r/2)} on the pair side, e^{iKm} on the composite side, all states
# vacuum-projected, measured windows used verbatim with a hard zero (cross
# block) or the disconnected skeleton (pair block) outside.
# =============================================================================

"""
    composite_pair_cross(H, Ω, χs, φ, d; seps, hops = 4, center = nothing) -> (sX, hX)

Measure the CROSS matrix elements between localized composite creators and the
two-particle basis on a finite chain: for the composite `χ̂_a` placed with its
LEFT EDGE at site `c` and a pair at `(c+m, c+m+r)`, return complex arrays
indexed `[a, ri, m + hops + 1]` with `r = seps[ri]`, `m ∈ -hops:hops`:

    sX = ⟨χ_a, c | (c+m, c+m+r)⟩ ,   hX = ⟨χ_a, c | (H − E₀) | (c+m, c+m+r)⟩ ,

where both bra and ket are vacuum-projected (`Pv = v − ⟨Ω,v⟩Ω`, exactly as in
[`pair_couplings`](@ref)). The composites in `χs` may have different supports.

Left-edge convention: pinning `χ̂_a` by its left edge (rather than its center)
is a pure bookkeeping choice — any fixed offset Δ only multiplies the
composite Bloch row by the overall phase `e^{iKΔ}`, i.e. rescales one basis
vector by a unimodular constant, which is irrelevant to the generalized
spectrum. What matters is that [`composite_pair_blocks`](@ref) consumes the
tables with the SAME convention, which it does.

For the composite–composite symbols do NOT remeasure anything here: call
`bulk_couplings(H, Ω, χs, d; rmax = hops)` — composites are one-body-like
channels, so the B₁ machinery applies verbatim and returns the
`N × N × (2 hops + 1)` arrays `S0, H0c` expected downstream.

`seps` defaults to `ℓ : ℓ+5` (`ℓ` = pair-creator support); all supports must
fit in `[1, L]` (`center` is chosen automatically mid-chain and can be
overridden).
"""
function composite_pair_cross(H::AbstractMatrix, Ω::AbstractVector, χs::Vector,
                              φ::AbstractMatrix, d::Integer;
                              seps = nothing, hops::Integer = 4,
                              center::Union{Nothing, Integer} = nothing)
    L = get_length_from_localdim(size(H, 1), d)
    ℓ = get_length_from_localdim(size(φ, 1), d)
    ℓχ = [get_length_from_localdim(size(χ, 1), d) for χ in χs]
    seps = seps === nothing ? collect(ℓ:(ℓ + 5)) : collect(seps)
    minimum(seps) >= ℓ || error("separations must be ≥ the creator support ℓ = $ℓ (disjoint supports)")
    M = hops
    smax = maximum(seps)
    ℓχmax = maximum(ℓχ)
    cmax = min(L - M - smax - ℓ + 1, L - ℓχmax + 1)
    c = center === nothing ? clamp((L - smax - ℓ + 2 - 2M) ÷ 2 + M, M + 1, max(M + 1, cmax)) :
                             Int(center)
    (c - M >= 1 && c + M + smax + ℓ - 1 <= L) ||
        error("pair window does not fit (need L ≥ $(2M + smax + ℓ)): L=$L, center=$c, smax=$smax, hops=$M, support=$ℓ")
    (c + ℓχmax - 1 <= L) ||
        error("composite support does not fit: left edge c=$c, largest composite support ℓχ=$ℓχmax, L=$L")
    E0 = real(dot(Ω, H * Ω))
    function pairket(j1, j2)
        v = embed_operator(φ, d, L, j2) * Ω
        v = embed_operator(φ, d, L, j1) * v
        v -= dot(Ω, v) * Ω
        v
    end
    nχ = length(χs)
    ns = length(seps)
    bras = map(1:nχ) do a
        v = embed_operator(χs[a], d, L, c) * Ω
        v -= dot(Ω, v) * Ω
        v
    end
    Hbras = [H * v - E0 * v for v in bras]
    sX = zeros(ComplexF64, nχ, ns, 2M + 1)
    hX = zeros(ComplexF64, nχ, ns, 2M + 1)
    for (rj, r) in enumerate(seps), m in -M:M
        v = pairket(c + m, c + m + r)
        for a in 1:nχ
            sX[a, rj, m + M + 1] = dot(bras[a], v)
            hX[a, rj, m + M + 1] = dot(Hbras[a], v)
        end
    end
    sX, hX
end

"""
    composite_pair_blocks(K, S0, H0c, sX, hX, s1, h1, S2, H2, seps;
                          hops = 4, Rmax = 300) -> (SK, HK, labels, rgrid)

Assemble the ENLARGED fixed-total-momentum-`K` generalized pair problem: an
`(nχ + n) × (nχ + n)` overlap/Hamiltonian pair over the basis
`{|K, χ_a⟩} ∪ {|K, r⟩, r ∈ rgrid}` with `n = length(rgrid)`. Inputs:

- `S0, H0c` — composite–composite bulk symbols from
  `bulk_couplings(H, Ω, χs, d; rmax = hops)`; Bloch sum
  `S_ab(K) = Σ_m e^{iKm} S0[a, b, m + hops + 1]` (composites are one-body-like,
  so the plain one-body phase applies);
- `sX, hX` — cross tables from [`composite_pair_cross`](@ref); the cross row
  carries the phase `e^{iK(m + r/2)}` matching the pair center-of-mass
  convention `|K, r⟩ = Σ_c e^{iK(c + r/2)}|c, c+r⟩`;
- `s1, h1, S2, H2, seps` — the one-/two-body data of [`pair_blocks`](@ref);
  the pair–pair block is delegated to `pair_blocks` verbatim (measured window
  + disconnected skeleton).

Outside the measured cross window (`r ∉ seps`) the cross element is set to a
HARD ZERO. This is justified — unlike for the pair–pair block, no skeleton is
needed: the cross overlap requires BOTH pair members to sit on the compact
composite support, so it decays exponentially in `m` and in `r` with no
disconnected (non-decaying) part to subtract.

Both outputs are Hermitized. `labels` is a `BitVector` of length `nχ + n`
with `true` marking composite indices (the first `nχ`), `rgrid` is the
relative-coordinate grid of the pair sector. Diagonalize with a Löwdin cutoff:
exact linear dependence of a composite on the pair basis shows up as a Gram
null direction and is dropped harmlessly (see [`friedrichs_lee`](@ref)).
"""
function composite_pair_blocks(K::Real, S0::Array{ComplexF64, 3}, H0c::Array{ComplexF64, 3},
                               sX::Array{ComplexF64, 3}, hX::Array{ComplexF64, 3},
                               s1::AbstractVector, h1::AbstractVector,
                               S2::Array{ComplexF64, 3}, H2::Array{ComplexF64, 3},
                               seps; hops::Integer = 4, Rmax::Integer = 300)
    seps = collect(seps)
    M = hops
    nχ = size(S0, 1)
    size(S0, 3) == 2M + 1 ||
        error("S0/H0c window $(size(S0, 3)) ≠ 2 hops + 1 = $(2M + 1): measure bulk_couplings with rmax = hops")
    size(sX, 3) == 2M + 1 ||
        error("sX/hX window $(size(sX, 3)) ≠ 2 hops + 1 = $(2M + 1): measure composite_pair_cross with the same hops")
    size(sX, 1) == nχ || error("sX has $(size(sX, 1)) composites but S0 has $nχ")
    size(sX, 2) == length(seps) || error("sX has $(size(sX, 2)) separations but seps has $(length(seps))")
    Spp, Hpp, rgrid = pair_blocks(K, s1, h1, S2, H2, seps; hops = M, Rmax = Rmax)
    n = length(rgrid)
    ntot = nχ + n
    SK = zeros(ComplexF64, ntot, ntot)
    HK = zeros(ComplexF64, ntot, ntot)
    # composite–composite block: one-body Bloch sum
    for a in 1:nχ, b in 1:nχ, m in -M:M
        ph = exp(im * K * m)
        SK[a, b] += ph * S0[a, b, m + M + 1]
        HK[a, b] += ph * H0c[a, b, m + M + 1]
    end
    # cross block: pair CM phase, hard zero outside the measured window
    sepidx = Dict(r => i for (i, r) in enumerate(seps))
    for (j, r) in enumerate(rgrid)
        haskey(sepidx, r) || continue
        rj = sepidx[r]
        for a in 1:nχ, m in -M:M
            ph = exp(im * K * (m + r / 2))
            SK[a, nχ + j] += ph * sX[a, rj, m + M + 1]
            HK[a, nχ + j] += ph * hX[a, rj, m + M + 1]
        end
    end
    SK[nχ+1:end, 1:nχ] .= adjoint(SK[1:nχ, nχ+1:end])
    HK[nχ+1:end, 1:nχ] .= adjoint(HK[1:nχ, nχ+1:end])
    SK[nχ+1:end, nχ+1:end] .= Spp
    HK[nχ+1:end, nχ+1:end] .= Hpp
    labels = BitVector([i <= nχ for i in 1:ntot])
    (SK + SK') / 2, (HK + HK') / 2, labels, rgrid
end

# Free-pair kinematic window from the one-body symbols ε(k) = h(k)/s(k):
# extrema of ε(K/2 + q) + ε(K/2 − q) over q. (Same few lines as inside
# two_particle_spectrum; duplicated on purpose to keep pair_bands.jl frozen.)
function _free_pair_window(K::Real, s1::AbstractVector, h1::AbstractVector)
    r1b = (length(s1) - 1) ÷ 2
    εk(k) = real(sum(exp(-im * k * r) * h1[r + r1b + 1] for r in -r1b:r1b) /
                 sum(exp(-im * k * r) * s1[r + r1b + 1] for r in -r1b:r1b))
    pairE = [εk(K / 2 + q) + εk(K / 2 - q) for q in range(-π, π; length = 601)]
    extrema(pairE)
end

"""
    friedrichs_lee(K, S0, H0c, sX, hX, s1, h1, S2, H2, seps;
                   hops = 4, Rmax = 300, gram_tol = 1e-6,
                   bound_frac = 0.2, bound_weight = 0.7, edge_pad = 1e-3)
        -> NamedTuple

Ab-initio Friedrichs–Lee reduction of the enlarged fixed-`K` problem of
[`composite_pair_blocks`](@ref). Procedure:

1. Löwdin-orthonormalize the PAIR subspace alone (cutoff `gram_tol` on its
   Gram spectrum) and diagonalize → discrete continuum levels `eps[a]` with
   orthonormal states `|a⟩`.
2. Project each composite out of the pair subspace USING THE S METRIC,
   then orthonormalize the residuals `χ̃_i` sequentially among themselves.
   A residual whose norm² falls below `gram_tol` times the largest raw
   composite norm² is a composite already representable by pairs: it is
   dropped with a warning and recorded in `dropped` (indices into the χ list;
   `kept` holds the complement).
3. The Friedrichs–Lee data: bare channel energies `Echi[i] = ⟨χ̃_i|H|χ̃_i⟩`
   (real vector over kept composites), couplings `g[i, a] = ⟨χ̃_i|H|a⟩`
   (complex `n_kept × n_cont` matrix), and the full channel block
   `Hcc[i, j] = ⟨χ̃_i|H|χ̃_j⟩` (equal to `Diagonal(Echi)` when a single
   composite is kept). Feed one row of `g` with `Echi[i]` and `eps` to
   [`resonance_parameters`](@ref) for the resonance position/width and the
   smooth resonant phase.
4. Independently, the FULL enlarged problem is diagonalized (global Löwdin
   cutoff `gram_tol`) → `energies`, with `isbound` from the same kinematic
   classifier as [`two_particle_spectrum`](@ref): the energy must lie outside
   the free-pair window `[lo(K) − edge_pad, hi(K) + edge_pad]` AND more than
   `bound_weight` of the (orthonormalized-frame) weight must be
   short-distance — on any composite or within the first `bound_frac`
   fraction of `rgrid` (composites are short-distance objects by
   construction).

Returns `(; Echi, g, Hcc, eps, dropped, kept, energies, isbound, rgrid)`.
"""
function friedrichs_lee(K::Real, S0::Array{ComplexF64, 3}, H0c::Array{ComplexF64, 3},
                        sX::Array{ComplexF64, 3}, hX::Array{ComplexF64, 3},
                        s1::AbstractVector, h1::AbstractVector,
                        S2::Array{ComplexF64, 3}, H2::Array{ComplexF64, 3},
                        seps; hops::Integer = 4, Rmax::Integer = 300,
                        gram_tol::Real = 1e-6, bound_frac::Real = 0.2,
                        bound_weight::Real = 0.7, edge_pad::Real = 1e-3)
    SK, HK, labels, rgrid = composite_pair_blocks(K, S0, H0c, sX, hX, s1, h1, S2, H2, seps;
                                                  hops = hops, Rmax = Rmax)
    nχ = count(labels)
    n = length(rgrid)
    p = (nχ + 1):(nχ + n)
    # --- 1. continuum: Löwdin + diagonalization of the pair subspace alone
    Spp = SK[p, p]
    Hpp = HK[p, p]
    eSp, USp = eigen(Hermitian(Spp))
    keepp = real.(eSp) .> gram_tol * maximum(real.(eSp))
    Wp = USp[:, keepp] * Diagonal(1 ./ sqrt.(real.(eSp[keepp])))
    εc, Vc = eigen(Hermitian(Wp' * Hpp * Wp))
    C = Wp * Vc                       # pair-grid coordinates of orthonormal |a⟩
    ncont = size(C, 2)
    # --- 2. residual composites (S-metric projection + sequential Gram–Schmidt)
    scale = maximum(real(SK[i, i]) for i in 1:nχ)   # raw composite norms² set the drop scale
    kept = Int[]
    dropped = Int[]
    R = Vector{ComplexF64}[]          # orthonormal residual coords in the enlarged basis
    for i in 1:nχ
        v = zeros(ComplexF64, nχ + n)
        v[i] = 1
        t = C' * SK[p, i]             # ⟨a|χ_i⟩ over the orthonormal continuum
        v[p] .-= C * t
        for w in R
            v .-= w .* dot(w, SK * v)
        end
        nrm2 = real(dot(v, SK * v))
        if nrm2 <= gram_tol * scale
            @warn "composite $i lies (numerically) inside the pair span at K = $K — dropped from the Friedrichs–Lee channel set" residual_norm² = nrm2 scale = scale
            push!(dropped, i)
        else
            v ./= sqrt(nrm2)
            push!(R, v)
            push!(kept, i)
        end
    end
    nk = length(R)
    # --- 3. Friedrichs–Lee data
    Echi = Float64[real(dot(w, HK * w)) for w in R]
    Hcc = ComplexF64[dot(R[i], HK * R[j]) for i in 1:nk, j in 1:nk]
    HC = HK[:, p] * C                 # H|a⟩ in enlarged coordinates
    g = ComplexF64[dot(R[i], HC[:, a]) for i in 1:nk, a in 1:ncont]
    # --- 4. full enlarged spectrum + kinematic bound classifier
    eS, US = eigen(Hermitian(SK))
    keepf = real.(eS) .> gram_tol * maximum(real.(eS))
    W = US[:, keepf] * Diagonal(1 ./ sqrt.(real.(eS[keepf])))
    E, V = eigen(Hermitian(W' * HK * W))
    lo, hi = _free_pair_window(K, s1, h1)
    ncut = max(1, round(Int, bound_frac * n))
    isbound = falses(length(E))
    for a in eachindex(E)
        (E[a] < lo - edge_pad || E[a] > hi + edge_pad) || continue
        w = W * V[:, a]
        w ./= norm(w)
        shortw = sum(abs2, @view w[1:nχ]) + sum(abs2, @view w[nχ .+ (1:ncut)])
        isbound[a] = shortw > bound_weight
    end
    (; Echi, g, Hcc, eps = real.(εc), dropped, kept,
       energies = real.(E), isbound, rgrid)
end

"""
    resonance_parameters(Echi, g, eps; eta = nothing, ngrid = 401) -> NamedTuple

Resonance position, width and smooth resonant phase from DISCRETE
Friedrichs–Lee data: a bare channel energy `Echi`, couplings `g[a]` to the
discrete continuum levels `eps[a]` (one kept composite of
[`friedrichs_lee`](@ref), or any single-channel FL model). The discrete sums
are smoothed on the scale `eta` (default: five mean level spacings of the
continuum window — a few spacings are needed for the kernel estimates to be
smooth, and the residual O(eta²) bias is negligible for dense spectra):

    Re Σ(E) = Σ_a |g_a|² (E − eps_a) / ((E − eps_a)² + eta²)         (smoothed principal value)
    Γ(E)    = 2π Σ_a |g_a|² N(E; eps_a, eta)                          (Gaussian-kernel estimate of
                                                                       2π ρ-weighted |g|², i.e. −2 Im Σ)

- `E_R` solves the resonance condition `E = Echi + Re Σ(E)` (bisection on the
  root nearest `Echi` inside the continuum window `[min eps, max eps]`, with a
  damped fixed-point fallback).
- If `Echi` lies OUTSIDE the continuum window the channel is a BOUND (or
  antibound-free) state: `Gamma = 0` is returned and `E_R` is the real root of
  the same condition, found by damped fixed-point iteration from `Echi` — the
  discrete analogue of the pole leaving the cut.
- `delta` is the smooth resonant phase sampled on `Egrid` (a uniform grid over
  the continuum window):  `delta(E) = atan2(Γ(E)/2, Echi + Re Σ(E) − E)`,
  rising through π/2 at the resonance. NOTE: this is the Friedrichs–Lee
  RESONANT part only — the potential-scattering (background) phase of the pair
  sector is not included.

Returns `(; E_R, Gamma, delta, Egrid, eta)` with `Gamma = Γ(E_R)` (or `0.0`
in the bound case).
"""
function resonance_parameters(Echi::Real, g::AbstractVector, eps::AbstractVector;
                              eta::Union{Nothing, Real} = nothing, ngrid::Integer = 401)
    isempty(eps) && error("empty continuum: eps must contain the discrete continuum levels")
    length(g) == length(eps) || error("g ($(length(g))) and eps ($(length(eps))) must have equal length")
    g2 = abs2.(g)
    emin, emax = extrema(eps)
    N = length(eps)
    η = eta === nothing ?
        max(5 * (emax - emin) / max(N - 1, 1), 1e-12 * max(1.0, abs(Echi))) : Float64(eta)
    ReΣ(E) = sum(g2[a] * (E - eps[a]) / ((E - eps[a])^2 + η^2) for a in 1:N)
    Γ(E) = 2π * sum(g2[a] * exp(-(E - eps[a])^2 / (2η^2)) for a in 1:N) / (sqrt(2π) * η)
    f(E) = Echi + ReΣ(E) - E
    fixedpoint(x0) = begin                       # damped fixed point, robust fallback
        E = float(x0)
        for _ in 1:300
            E = 0.5 * E + 0.5 * (Echi + ReΣ(E))
        end
        E
    end
    local E_R::Float64, Gamma::Float64
    if emin <= Echi <= emax
        # resonance: bisect on the sign change nearest to Echi inside the window
        Es = range(emin, emax; length = 2001)
        fs = [f(E) for E in Es]
        E_R = NaN
        best = Inf
        for i in 1:(length(Es) - 1)
            sign(fs[i]) != sign(fs[i + 1]) || continue
            a, b = Es[i], Es[i + 1]
            fa = fs[i]
            for _ in 1:60
                mid = (a + b) / 2
                if sign(f(mid)) == sign(fa)
                    a = mid
                else
                    b = mid
                end
            end
            cand = (a + b) / 2
            if abs(cand - Echi) < best
                best = abs(cand - Echi)
                E_R = cand
            end
        end
        isnan(E_R) && (E_R = fixedpoint(Echi))
        Gamma = Γ(E_R)
    else
        # bound state: the pole sits on the real axis outside the continuum
        E_R = fixedpoint(Echi)
        Gamma = 0.0
    end
    Egrid = collect(range(emin, emax; length = ngrid))
    delta = [atan(Γ(E) / 2, Echi + ReΣ(E) - E) for E in Egrid]
    (; E_R, Gamma, delta, Egrid, eta = η)
end
