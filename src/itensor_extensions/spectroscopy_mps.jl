# =============================================================================
# spectroscopy_mps.jl — MPS port of the bulk-reconstruction spectroscopy
#
# The OBC production mode of the method: everything the matrix pipeline does
# (trap seeds → transition-operator creators → bulk couplings → Wannier
# bands) expressed in MPS/MPO operations, so it scales to L ~ 24–30 and
# beyond. The Fourier/linear-algebra stage (`wannier_bands`) is shared with
# the matrix backend unchanged.
#
# WINDOW CONVENTION (shared with the matrix backend): a dense window operator
# A ∈ C^{d^ℓ × d^ℓ} on MPS sites w1 < w2 < … < wℓ indexes its rows/columns in
# kron order with the FIRST window site slowest (i.e. row = i_{w1} ⊗ … ⊗
# i_{wℓ}, last site fastest) — the same convention as `transition_operator`
# on the matrix side. All conversions below respect it.
# =============================================================================

# window ITensor with the shared ordering: rows (output) primed.
function _window_gate(A::AbstractMatrix, sw::Vector{<:Index})
    ℓ = length(sw)
    d = dim(sw[1])
    size(A) == (d^ℓ, d^ℓ) || error("window operator size $(size(A)) ≠ ($(d^ℓ), $(d^ℓ))")
    T = reshape(Matrix{ComplexF64}(A), ntuple(_ -> d, 2ℓ))
    # Julia reshape is column-major: dims are (r_wℓ … r_w1, c_wℓ … c_w1),
    # fastest first — so indices are listed in REVERSED window order.
    ITensor(T, (prime.(reverse(sw))..., reverse(sw)...))
end

"""
    mps_transition_operator(ψ::MPS, Ω::MPS, window::UnitRange) -> Matrix{ComplexF64}

The reduced transition operator `A = Tr_rest |ψ⟩⟨Ω|` on the contiguous
`window` of sites, computed by exact transfer-environment contraction
(left environment up to the window, right environment after it, window
physical indices left open on both ket and bra).

Follows the shared window convention (first window site slowest), so the
result is interchangeable with the matrix-side [`transition_operator`](@ref)
and feeds [`mps_apply_window`](@ref) / [`mps_bulk_couplings`](@ref) directly.
`ψ` and `Ω` may carry different site indices (they are matched positionally).
Cost: O(L · χψ χΩ d) plus the open-window blowup d^{2ℓ} — fine for ℓ ≲ 6.
"""
function mps_transition_operator(ψ::MPS, Ω::MPS, window::UnitRange{<:Integer})
    L = length(ψ)
    length(Ω) == L || error("ψ and Ω must have the same length")
    1 <= first(window) <= last(window) <= L || error("window $window outside 1:$L")
    Ωc = copy(Ω)
    substitute_siteinds!(Ωc, ψ)
    s = siteinds(ψ)
    E = ITensor(1.0)
    for n in 1:first(window)-1
        E = E * ψ[n] * dag(Ωc[n])
    end
    for n in window
        E = E * ψ[n] * prime(dag(Ωc[n]), s[n])
    end
    R = ITensor(1.0)
    for n in L:-1:last(window)+1
        R = R * ψ[n] * dag(Ωc[n])
    end
    T = E * R
    sw = [s[n] for n in window]
    d = dim(sw[1]); ℓ = length(sw)
    arr = Array(T, (reverse(sw)..., reverse(prime.(sw))...))
    Matrix{ComplexF64}(reshape(arr, d^ℓ, d^ℓ))
end

"""
    mps_apply_window(A::AbstractMatrix, Ω::MPS, site::Integer;
                     cutoff = default_cutoff, maxdim = default_maxdim) -> MPS

Apply a dense `ℓ`-site window operator `A` (shared window convention) at
position `site` (leftmost window site) of the MPS `Ω`: build the window MPO,
pad it with identities ([`insert_local`](@ref)), match the site indices and
`apply`. Returns the (unnormalized) resulting MPS.
"""
function mps_apply_window(A::AbstractMatrix, Ω::MPS, site::Integer;
                          cutoff = default_cutoff, maxdim = default_maxdim)
    L = length(Ω)
    d = dim(siteinds(Ω)[1])
    ℓ = round(Int, log(d, size(A, 1)))
    d^ℓ == size(A, 1) || error("window operator size is not a power of the local dimension")
    1 <= site <= L - ℓ + 1 || error("window of $ℓ sites does not fit at site $site of $L")
    sw = siteinds(d, ℓ)
    win = MPO(_window_gate(A, sw), sw)
    full = insert_local(site - 1, win, L - site - ℓ + 1)
    replace_siteinds!(full, siteinds(Ω))
    apply(full, Ω; cutoff, maxdim)
end

"""
    deformed_hamiltonian_mpo(H0::AbstractMatrix, d, L;
                             weight = parabolic_weight(L),
                             cutoff = default_cutoff, maxdim = default_maxdim) -> MPO

MPO of the deformed (trap) Hamiltonian `H' = Σ_j weight(j) h_j` on an OPEN
chain of `L` sites, from the local block `H0` — the MPS counterpart of the
matrix [`deformed_hamiltonian`](@ref) with `pbc = false`. `weight(j)` is the
coefficient of the block starting at site `j` (1-based), defaulting to the
open-chain parabolic trap. Built by weighted summation of translated local
MPOs ([`summation_local`](@ref) with `convolution = weight`).
"""
function deformed_hamiltonian_mpo(H0::AbstractMatrix, d::Integer, L::Integer;
                                  weight = parabolic_weight(L),
                                  cutoff = default_cutoff, maxdim = default_maxdim)
    h0 = mpo_from_matrix(Matrix(H0), d; cutoff, maxdim)
    summation_local(h0, L; convolution = weight, cutoff, maxdim)
end

"""
    dmrg_seeds(Hmpo::MPO, sites, N;
               nsweeps = 12, maxdim = [10, 20, 40, 80, 120],
               cutoff = 1e-10, penalty = 20.0, linkdims = 8,
               outputlevel = 0) -> (energies, Vector{MPS})

The `N` lowest eigenstates of an MPO (typically the trap Hamiltonian) by
successive DMRG runs with an orthogonality penalty against the previously
found states — the MPS counterpart of [`trap_seeds`](@ref). Returns the
energies and the states, lowest first.

DMRG starting states are drawn from the global RNG — seed it
(`Random.seed!`) for reproducible seeds when the trap spectrum is
near-degenerate.
"""
function dmrg_seeds(Hmpo::MPO, sites, N::Integer;
                    nsweeps::Integer = 12, maxdim = [10, 20, 40, 80, 120],
                    cutoff::Real = 1e-10, penalty::Real = 20.0,
                    linkdims::Integer = 8, outputlevel::Integer = 0)
    energies = Float64[]; states = MPS[]
    for n in 1:N
        ψ0 = random_mps(sites; linkdims)
        E, ψ = if isempty(states)
            dmrg(Hmpo, ψ0; nsweeps, maxdim, cutoff, outputlevel)
        else
            dmrg(Hmpo, states, ψ0; nsweeps, maxdim, cutoff, outputlevel, weight = penalty)
        end
        push!(energies, E); push!(states, ψ)
    end
    energies, states
end

"""
    mps_bulk_couplings(Hmpo::MPO, Ω::MPS, creators, center, rmax;
                       cutoff = default_cutoff, maxdim = default_maxdim)
        -> (S, Hc)

Bulk matrix elements between translated dressed states
`|α, j⟩ = P_Ω⊥ φ̂_α^j |Ω⟩` on an MPS — the production (OBC / tensor-network)
counterpart of the matrix [`bulk_couplings`](@ref), with identical output
format: `N × N × (2 rmax + 1)` arrays indexed `[α, β, r + rmax + 1]`,

    S[α,β,·]  = ⟨α, center | β, center + r⟩ ,
    Hc[α,β,·] = ⟨α, center | (H − E0) | β, center + r⟩ ,

`E0 = ⟨Ω|H|Ω⟩`. The vacuum projection is performed ALGEBRAICALLY (never
building projected MPS): with `u = φ^c Ω`, `v = φ^{c+r} Ω`, `a = ⟨Ω|u⟩`,
`b = ⟨Ω|v⟩`,

    ⟨ũ|ṽ⟩          = ⟨u|v⟩ − ā b ,
    ⟨ũ|H − E0|ṽ⟩   = ⟨u|H|v⟩ − E0⟨u|v⟩ − b (⟨u|H|Ω⟩ − E0 ā)
                                       − ā (⟨Ω|H|v⟩ − E0 b) .

Feed the output straight into [`wannier_bands`](@ref). Choose `center` and
`rmax` so all windows stay in the translation-invariant bulk.
"""
function mps_bulk_couplings(Hmpo::MPO, Ω::MPS, creators::Vector{<:AbstractMatrix},
                            center::Integer, rmax::Integer;
                            cutoff = default_cutoff, maxdim = default_maxdim)
    N = length(creators)
    L = length(Ω)
    d = dim(siteinds(Ω)[1])
    ℓs = [round(Int, log(d, size(A, 1))) for A in creators]
    (center - rmax >= 1 && center + rmax + maximum(ℓs) - 1 <= L) ||
        error("window [center ± rmax + support] exceeds the chain; reduce rmax or enlarge L")
    E0 = real(inner(Ω', Hmpo, Ω))
    positions = (center - rmax):(center + rmax)
    kets = Dict{Tuple{Int, Int}, MPS}()
    va = Dict{Tuple{Int, Int}, ComplexF64}()      # ⟨Ω|ket⟩
    vh = Dict{Tuple{Int, Int}, ComplexF64}()      # ⟨Ω|H|ket⟩
    for α in 1:N, p in positions
        k = mps_apply_window(creators[α], Ω, p; cutoff, maxdim)
        kets[(α, p)] = k
        va[(α, p)] = inner(Ω, k)
        vh[(α, p)] = inner(Ω', Hmpo, k)
    end
    S = zeros(ComplexF64, N, N, 2rmax + 1)
    Hc = zeros(ComplexF64, N, N, 2rmax + 1)
    for β in 1:N, r in -rmax:rmax
        v = kets[(β, center + r)]
        b = va[(β, center + r)]
        hv = vh[(β, center + r)]                  # ⟨Ω|H|v⟩
        for α in 1:N
            u = kets[(α, center)]
            a = va[(α, center)]
            hu = vh[(α, center)]                  # ⟨Ω|H|u⟩ ⇒ ⟨u|H|Ω⟩ = conj(hu)
            Suv = inner(u, v)
            Huv = inner(u', Hmpo, v)
            S[α, β, r + rmax + 1] = Suv - conj(a) * b
            Hc[α, β, r + rmax + 1] = Huv - E0 * Suv -
                b * (conj(hu) - E0 * conj(a)) -
                conj(a) * (hv - E0 * b)
        end
    end
    S, Hc
end
