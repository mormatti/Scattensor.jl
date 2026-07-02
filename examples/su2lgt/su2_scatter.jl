# =============================================================================
# su2_scatter.jl  --  baryon dispersion + baryon-baryon Lüscher phase shift
# for the periodic (1+1)D SU(2) lattice gauge theory (arXiv:2511.00154).
#
# The lightest particle in the Fig-4 strong-coupling regime is the BARYON
# (color-singlet diquark, B=1, ~2m, no confining string).  We extract:
#   * the single-baryon dispersion ε(k)            from the B = 1 sector,
#   * the two-baryon finite-volume levels          from the B = 2 sector,
# and feed them into the same Lüscher machinery used for the spin chains.
#
# Baryon number  B = ¼ Σ σ^z  ⇒  #up = N + 2B.   Translation = 4-qubit (2-site)
# shift  ⇒  N/2 momentum cells, K = 2π·Kint/(N/2).
# =============================================================================

const HERE = @__DIR__
include(joinpath(HERE, "su2_lgt.jl"))
using LinearAlgebra

# basis indices (1-based) of the B-sector:  Σσ^z = 4B  ⇔  #up = N + 2B
function bsector_indices(N::Int, B::Int)
    nq = 2N; want = N + 2B
    [s+1 for s in 0:(2^nq-1) if count_ones(s) == want]
end

# 4-qubit (2-site) cyclic translation as a permutation on 2^nq basis states.
# Qubit j is bit (nq-j) of the state index (qubit 1 = most significant, per op_at).
function translation_4(nq::Int)
    perm = Vector{Int}(undef, 2^nq)
    for s in 0:(2^nq-1)
        t = 0
        for j in 1:nq
            bit = (s >> (nq-j)) & 1           # value of qubit j
            jt = mod(j-1+4, nq) + 1           # qubit j → qubit j+4
            t |= bit << (nq-jt)
        end
        perm[s+1] = t+1
    end
    return sparse(1:2^nq, perm, ones(2^nq))   # T|s⟩ = |perm(s)⟩
end

# energy levels labelled by lattice momentum within a sector.
# Robust: project onto each momentum sector with the Hermitian projector
#   P_K = (1/ncells) Σ_r e^{-iθr} T^r ,  θ = 2π Kint/ncells,
# orthonormalise its range, and diagonalise H inside it.
# Returns Vector of (Kint, K, energy), energy-sorted per Kint.
function sector_momentum_levels(H, T4, idx::Vector{Int}, ncells::Int)
    Hs = Matrix(H[idx, idx]); Ts = Matrix(T4[idx, idx])
    @assert norm(Hs - Hs') < 1e-8 "non-Hermitian H"
    @assert norm(Ts*Ts' - I) < 1e-8 "non-unitary T"
    @assert norm(Hs*Ts - Ts*Hs) < 1e-6 "[H,T]≠0"
    Tpow = [Ts^r for r in 0:ncells-1]
    out = NamedTuple[]
    for Kint in 0:ncells-1
        θ = 2π*Kint/ncells
        P = sum(cis(-θ*r) .* Tpow[r+1] for r in 0:ncells-1) ./ ncells
        P = (P + P')/2
        ep = eigen(Hermitian(P))
        cols = findall(>(0.5), real.(ep.values))   # eigenvalue ≈ 1 ⇒ in this sector
        isempty(cols) && continue
        Q = ep.vectors[:, cols]                     # orthonormal basis of the K-block
        HK = Hermitian(Q' * Hs * Q)
        for e in eigvals(HK)
            push!(out, (Kint=Kint, K=2π*Kint/ncells, energy=e))
        end
    end
    sort!(out, by = r -> (r.Kint, r.energy))
    return out
end

# ---- validate: single-baryon dispersion (B=1) ------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    for N in (6,)
        H, nq = build_su2_hamiltonian(N; a=1.0, g=5.0, m=0.2)
        T4 = translation_4(nq)
        ncells = N ÷ 2
        # vacuum energy from B=0
        lev0 = sector_momentum_levels(H, T4, bsector_indices(N,0), ncells)
        E0 = minimum(r.energy for r in lev0)
        # B=1 single-baryon band: lowest level in each momentum sector
        lev1 = sector_momentum_levels(H, T4, bsector_indices(N,1), ncells)
        println("\nN=$N  (ncells=$ncells, E0=$(round(E0,digits=4)))")
        println("  single-baryon band ε(K)=E_min(B=1,K) − E0:")
        for Kint in 0:ncells-1
            es = sort([r.energy for r in lev1 if r.Kint==Kint])
            isempty(es) && continue
            println("    Kint=$Kint  K=$(round(2π*Kint/ncells,digits=3))  ε=$(round(es[1]-E0,digits=4))")
        end
    end
end
