# =============================================================================
# sectors.jl — symmetry-sector decomposition of a Hamiltonian (exact ED speedup)
#
# Given a Hamiltonian H and a MONOMIAL unitary symmetry S (one nonzero entry
# per column — permutation matrices possibly dressed with phases; lattice
# translations, spin flips, particle-hole, etc.), build the isometries onto
# each symmetry eigenspace from the ORBITS of the underlying permutation.
# The projected blocks have dimension ~ dim(H)/order instead of dim(H):
# diagonalizing the blocks is dramatically cheaper than projecting in the
# full space (the strategy of `dispersion_relation`), with identical spectra.
#
# Main entry points:
#   symmetry_blocks(S; order)          → sector isometries for any monomial S
#   momentum_basis(T, L)               → the translation case, labeled by k/π
#   sector_spectrum(H, T, L; nlevels)  → fast drop-in for dispersion_relation
# =============================================================================

"""
    SymmetrySector

One symmetry eigenspace of a monomial unitary `S`.

# Fields
- `label::Rational{Int}`: the eigenvalue of `S` is `exp(iπ · label)`. For a
  lattice translation this is `k/π`; for a Z₂ symmetry it is `0//1` or `1//1`.
- `basis::SparseMatrixCSC{ComplexF64,Int}`: isometry `B` (dim(H) × sector dim)
  whose columns are the symmetry-adapted (Bloch) states: `S B = exp(iπ·label) B`
  and `B† B = I`.
"""
struct SymmetrySector
    label::Rational{Int}
    basis::SparseMatrixCSC{ComplexF64, Int}
end

# --- monomial structure extraction ------------------------------------------
# S must have exactly one nonzero per column: S e_i = w[i] e_{p[i]}.
function _monomial_structure(S::AbstractMatrix)
    n = size(S, 1)
    size(S, 2) == n || error("symmetry operator must be square")
    Ssp = sparse(S)
    p = zeros(Int, n); w = zeros(ComplexF64, n)
    rows = rowvals(Ssp); vals = nonzeros(Ssp)
    for i in 1:n
        rng = nzrange(Ssp, i)
        idx = [r for r in rng if abs(vals[r]) > 1e-13]
        length(idx) == 1 || error("operator is not monomial: column $i has $(length(idx)) nonzeros")
        p[i] = rows[idx[1]]
        w[i] = vals[idx[1]]
        abs(abs(w[i]) - 1) < 1e-10 || error("operator is not unitary-monomial: |entry| = $(abs(w[i])) in column $i")
    end
    p, w
end

"""
    symmetry_blocks(S::AbstractMatrix; order::Integer) -> Vector{SymmetrySector}

Decompose the Hilbert space into eigenspaces of a **monomial** unitary `S`
with `S^order ∝ I` (lattice translation: `order = L`; spin flip: `order = 2`).

The construction is orbit-based: the permutation underlying `S` is split into
cycles; each cycle of length `ℓ` with accumulated phase `W` (product of the
monomial phases around the cycle) contributes one symmetry-adapted vector for
every `ℓ`-th root `χ` of `W`,

    v_χ = ℓ^{-1/2} Σ_{r=0}^{ℓ-1} χ^{-r} W_r |p^r(i)⟩ ,     S v_χ = χ v_χ .

Sectors are labeled by `label = arg(χ)/π` rounded to the rational grid
`2m//order` (translations: `k/π`). Only `label ∈ [0, 1]` sectors are returned
when `S` is real (the `-k` partners are complex conjugates); pass
`fullzone = true` to keep all.

Returns the sectors sorted by label. The direct sum of all sector dimensions
equals `dim(H)` (with `fullzone = true`).
"""
function symmetry_blocks(S::AbstractMatrix; order::Integer, fullzone::Bool = false)
    n = size(S, 1)
    p, w = _monomial_structure(S)
    visited = falses(n)
    # triplet buffers per label index m (χ = exp(2πi m / order))
    rowbuf = [Int[] for _ in 1:order]
    colbuf = [Int[] for _ in 1:order]
    valbuf = [ComplexF64[] for _ in 1:order]
    ncols = zeros(Int, order)
    for i0 in 1:n
        visited[i0] && continue
        # walk the cycle, accumulating phases
        cyc = Int[]; phs = ComplexF64[]
        i = i0; W = 1.0 + 0im
        while true
            push!(cyc, i); push!(phs, W)
            visited[i] = true
            W *= w[i]
            i = p[i]
            i == i0 && break
        end
        ℓ = length(cyc)
        mod(order, ℓ) == 0 || error("cycle length $ℓ does not divide the symmetry order $order")
        # allowed eigenvalues: χ^ℓ = W (phase closure)
        θW = angle(W)
        for j in 0:ℓ-1
            χ = exp(im * (θW + 2π * j) / ℓ)
            # map χ onto the global character grid exp(2πi m / order)
            m = mod(round(Int, angle(χ) * order / (2π)), order)
            ncols[m + 1] += 1
            c = ncols[m + 1]
            for r in 1:ℓ
                push!(rowbuf[m + 1], cyc[r])
                push!(colbuf[m + 1], c)
                push!(valbuf[m + 1], conj(χ)^(r - 1) * phs[r] / sqrt(ℓ))
            end
        end
    end
    sectors = SymmetrySector[]
    for m in 0:order-1
        ncols[m + 1] == 0 && continue
        f = mod(2m // order + 1, 2) - 1                # label in (-1, 1]
        (!fullzone && f < 0) && continue
        B = sparse(rowbuf[m + 1], colbuf[m + 1], valbuf[m + 1], n, ncols[m + 1])
        push!(sectors, SymmetrySector(f, B))
    end
    sort!(sectors, by = s -> s.label)
    sectors
end

"""
    momentum_basis(T::AbstractMatrix, L::Integer; fullzone=false) -> Vector{SymmetrySector}

Momentum-sector isometries from the translation operator (`order = L`).
Labels are `k/π ∈ [0, 1]` (or the full zone if `fullzone = true`).
"""
momentum_basis(T::AbstractMatrix, L::Integer; fullzone::Bool = false) =
    symmetry_blocks(T; order = L, fullzone)

"""
    sector_hamiltonian(H::AbstractMatrix, sec::SymmetrySector) -> SparseMatrixCSC

The Hamiltonian block `B† H B` in the symmetry sector (Hermitized).
"""
function sector_hamiltonian(H::AbstractMatrix, sec::SymmetrySector)
    B = sec.basis
    Hk = B' * sparse(H) * B
    (Hk + Hk') / 2
end

"""
    sector_spectrum(H, T, L; nlevels=2, dense_cutoff=4096, fullzone=false,
                    lift=true, verbose=false) -> Vector{BlochState}

Fast, exact replacement for [`dispersion_relation`](@ref): decompose into
momentum sectors of dimension ~`dim(H)/L` via [`momentum_basis`](@ref),
diagonalize each block (densely below `dense_cutoff`, Krylov above), and
return the `nlevels` lowest [`BlochState`](@ref)s per sector.

With `lift = true` (default) the eigenvectors are lifted back to the full
Hilbert space (`B v`), so downstream code (overlaps, localizability, …) works
unchanged. Set `lift = false` to keep sector-internal vectors and save memory.
"""
function sector_spectrum(H::AbstractMatrix, T::AbstractMatrix, L::Integer;
                         nlevels::Integer = 2, dense_cutoff::Integer = 4096,
                         fullzone::Bool = false, lift::Bool = true,
                         verbose::Bool = false)
    states = BlochState[]
    for sec in momentum_basis(T, L; fullzone)
        Hk = sector_hamiltonian(H, sec)
        nk = size(Hk, 1)
        verbose && println("sector k/π = $(sec.label): dim $nk")
        en, st = if nk <= dense_cutoff
            e, V = eigen(Hermitian(Matrix(Hk)))
            e, collect(eachcol(V))
        else
            e, v, _ = eigsolve(Hk, nk, min(nlevels, nk), :SR, ComplexF64;
                               ishermitian = true,
                               krylovdim = max(2nlevels, nlevels + 20))
            real.(e), v
        end
        for i in 1:min(nlevels, length(en))
            ψ = lift ? sec.basis * st[i] : st[i]
            push!(states, BlochState(ψ, real(en[i]), sec.label))
        end
    end
    states
end
