# =============================================================================
# trap.jl — variational quasiparticles from deformed Hamiltonians
#
# Deform the Hamiltonian DENSITY (not the potential): H' = Σ_j V(j) h_j with a
# smooth PBC weight. The default V(j) = 1 − cos(2π(j−i_c)/L) is the sine-square
# deformation (SSD); every quasiparticle species self-traps at the weight
# minimum, so the lowest eigenstates of H' are localized packets of ALL low
# bands at once (they are the maximally-localized Wannier functions of each
# band — see the Spectroscopy notes for the Wannier/SSD/entanglement-
# Hamiltonian genealogy).
#
# Pipeline:  deformed_hamiltonian → (Krylov: N seeds) → effective_pair
#            → dispersion_relation(H̃, T̃, L)   [all low bands, one shot]
# =============================================================================

"""
    ssd_weight(L; center = L ÷ 2) -> Function

The sine-square-deformation weight `V(j) = 1 − cos(2π (j − center)/L)`:
smooth on the ring, zero (quadratic) at `center`, maximal at the antipode.
Natural for PBC traps.
"""
ssd_weight(L::Integer; center::Integer = L ÷ 2) = j -> 1.0 - cos(2π * (j - center) / L)

"""
    parabolic_weight(L; center = (L + 1) ÷ 2, amp = 2.0) -> Function

The open-chain trap weight `V(j) = amp · ((j − center)/(L/2))²`: zero
(quadratic) at `center`, growing toward the walls. Natural for OBC traps
(the PBC-smoothness of [`ssd_weight`](@ref) is not needed without a ring).
"""
parabolic_weight(L::Integer; center::Integer = (L + 1) ÷ 2, amp::Real = 2.0) =
    j -> amp * ((j - center) / (L / 2))^2

"""
    deformed_hamiltonian(H0, d, L; weight = ssd_weight(L), pbc = true, T = nothing)
        -> SparseMatrixCSC

Build `H' = Σ_j weight(j) · h_j` from the local Hamiltonian block `H0`
(acting on `L0 = log_d(size(H0,1))` sites), with `h_j` starting at site `j`
(1-based). This is the "trap": a modulation of the local Hamiltonian
density, NOT an external potential.

With `pbc = true` the blocks wrap around the ring (translation-conjugated
placement); with `pbc = false` only the `L − L0 + 1` fitting positions are
used — combine with [`parabolic_weight`](@ref) for the open-chain trap.
`T` can be passed to avoid rebuilding the translation operator (PBC only).
The result is Hermitized against accumulation of numerical noise.
"""
function deformed_hamiltonian(H0::AbstractMatrix, d::Integer, L::Integer;
                              weight = ssd_weight(L), pbc::Bool = true,
                              T::Union{Nothing, AbstractMatrix} = nothing)
    L0 = get_length_from_localdim(size(H0, 1), d)
    Hprime = spzeros(ComplexF64, d^L, d^L)
    for j in 1:(L - L0 + 1)
        Hprime += weight(j) * embed_operator(H0, d, L, j)
    end
    if pbc
        Top = T === nothing ? operator_translation(SparseMatrixCSC, d, L) : T
        Hcur = SparseMatrixCSC{ComplexF64}(embed_operator(H0, d, L, L - L0 + 1))
        for j in (L - L0 + 2):L
            Hcur = Top * Hcur * Top'                 # wrap terms across the boundary
            Hprime += weight(j) * Hcur
        end
    end
    (Hprime + Hprime') / 2
end

"""
    trap_seeds(Hprime, N; krylovdim = max(2N, N + 20)) -> (energies, states)

The `N` lowest eigenpairs of the deformed Hamiltonian (Krylov). The states
are localized packets of the low quasiparticle species, pinned at the trap
center — the variational seeds of the effective-band construction.
"""
function trap_seeds(Hprime::AbstractMatrix, N::Integer; krylovdim::Integer = max(2N, N + 20))
    en, st, _ = eigsolve(Hprime, size(Hprime, 1), N, :SR, ComplexF64;
                         ishermitian = true, krylovdim)
    real.(en[1:N]), st[1:N]
end

"""
    effective_pair(H, T, seeds, L; tol = 1e-8) -> (H̃, T̃, W)

From `N` localized seed states, build the translated (non-orthogonal) basis
`B = {T^p |n⟩ : n = 1..N, p = 0..L-1}`, Löwdin-orthogonalize with a spectral
cutoff on the overlap matrix (keep `S`-eigenvalues `> tol · max`), and return
the effective pair

    H̃ = W† (B†HB) W ,   T̃ = W† (B†TB) W ,   W = U_S Σ_S^{-1/2} ,

on which [`dispersion_relation`](@ref) reconstructs all the bands spanned by
the seeds, at cost `(N·L)³` instead of `d^L`. `W` maps effective vectors back
to basis-B coefficients (full-space lift: `Vmat · W · v`).
"""
function effective_pair(H::AbstractMatrix, T::AbstractMatrix, seeds::Vector, L::Integer;
                        tol::Real = 1e-8)
    dimH = size(H, 1)
    nL = Int(L)
    N = length(seeds)
    Vmat = zeros(ComplexF64, dimH, N * nL)
    for n in 1:N
        v = Vector{ComplexF64}(seeds[n])
        for p in 0:(nL - 1)
            Vmat[:, (n - 1) * nL + p + 1] = v
            v = T * v
        end
    end
    S = Vmat' * Vmat;        S = (S + S') / 2
    HB = Vmat' * (H * Vmat); HB = (HB + HB') / 2
    TB = Vmat' * (T * Vmat)
    eS, US = eigen(Hermitian(S))
    keep = real.(eS) .> tol * maximum(real.(eS))
    W = US[:, keep] * Diagonal(1 ./ sqrt.(real.(eS[keep])))
    H̃ = W' * HB * W; H̃ = (H̃ + H̃') / 2
    T̃ = W' * TB * W
    Matrix(H̃), Matrix(T̃), W
end
