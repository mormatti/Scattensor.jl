# =============================================================================
# su2_lgt.jl  --  (1+1)D SU(2) lattice gauge theory on a PERIODIC ring
# =============================================================================
#
# Spin/qubit Hamiltonian of arXiv:2511.00154 (gaugeless formulation, eq. 2.8),
# rewritten with PERIODIC boundary conditions so that lattice momentum is a good
# quantum number (needed for the Lüscher phase-shift extraction).
#
#   * N color-doubled sites n=1..N  →  2N qubits j=1..2N.
#     site n carries the two colors at qubits (2n-1, 2n):  φ_n = (ψ_{2n-1}, ψ_{2n}).
#   * Kinetic  :  (1/2a) Σ_j (σ⁺_j σ⁻_{j+1} + h.c.)          (n.n. hop, periodic)
#   * Mass     :  m a Σ_n ((-1)^n/2)(σ^z_{2n-1}+σ^z_{2n})    (staggered)
#   * Electric :  (a g²/2) Σ_{n,m} G(n-m) ( Q_n·Q_m )         (PERIODIC 1D Coulomb)
#         color charges   Q_n^z = ¼(σ^z_{2n-1}-σ^z_{2n}),
#                         Q_n^+ = σ⁺_{2n-1}σ⁻_{2n},  Q_n^- = (Q_n^+)† ,
#                         Q_n·Q_m = Q_n^z Q_m^z + ½(Q_n^+Q_m^- + Q_n^-Q_m^+),
#         G = (-Δ)^{-1}  the periodic lattice Green's function (zero mode removed),
#         which is the translation-invariant analogue of the paper's open-boundary
#         Σ_{k≤n} (linear confining) Coulomb potential.
#
# Conserved baryon number  B = ¼ Σ_j σ^z_j ;  we work in the B = 0 sector.
# Parameters (paper's Fig. 4 regime):  a = 1,  g a = 5,  m a = 0.2  (x=1/(ga)²=0.04).
# =============================================================================

using SparseArrays, LinearAlgebra

# ---- single-qubit Paulis on 2N qubits (qubit 1 = most significant) ----------
const _I2 = sparse([1.0 0; 0 1.0])
const _SZ = sparse([1.0 0; 0 -1.0])
const _SP = sparse([0.0 1.0; 0 0])    # σ⁺ = |↑⟩⟨↓|  (raising; |↑⟩=occupied)
const _SM = sparse([0.0 0; 1.0 0])    # σ⁻

# operator `op` acting on qubit j of nq qubits
function op_at(op, j::Int, nq::Int)
    m = sparse(1.0I, 1, 1)
    for q in 1:nq
        m = kron(m, q == j ? op : _I2)
    end
    return m
end

# periodic 1D lattice Green's function G(d) = (-Δ)^{-1}, zero mode removed
function periodic_green(N::Int)
    G = zeros(N)
    for d in 0:N-1
        s = 0.0
        for jj in 1:N-1
            q = 2π*jj/N
            s += cos(q*d) / (2 - 2cos(q))
        end
        G[d+1] = s/N
    end
    return G                          # G[d+1] = G(d), even & periodic
end

"""
    build_su2_hamiltonian(N; a, g, m) -> (H, nq)

Full 2^(2N) sparse Hamiltonian of the periodic SU(2) LGT.
"""
function build_su2_hamiltonian(N::Int; a=1.0, g=5.0, m=0.2, Λ=100.0)
    nq = 2N
    SZ(j)=op_at(_SZ,j,nq); SP(j)=op_at(_SP,j,nq); SM(j)=op_at(_SM,j,nq)
    H = spzeros(ComplexF64, 2^nq, 2^nq)

    # kinetic: (1/2a) Σ_j (σ⁺_j σ⁻_{j+1} + h.c.), periodic
    for j in 1:nq
        jp = j == nq ? 1 : j+1
        H += (1/(2a)) * (SP(j)*SM(jp) + SM(j)*SP(jp))
    end

    # staggered mass: m a Σ_n ((-1)^n/2)(σ^z_{2n-1}+σ^z_{2n})
    for n in 1:N
        H += (m*a) * ((-1)^n/2) * (SZ(2n-1) + SZ(2n))
    end

    # color charges per site
    Qz = [0.25*(SZ(2n-1) - SZ(2n)) for n in 1:N]
    Qp = [SP(2n-1)*SM(2n) for n in 1:N]
    Qm = [SM(2n-1)*SP(2n) for n in 1:N]

    # electric energy: (a g²/2) Σ_{n,m} G(n-m) (Q_n·Q_m)
    G = periodic_green(N)
    pref = a*g^2/2
    for n in 1:N, mm in 1:N
        d = mod(n-mm, N)
        gg = G[d+1]
        gg == 0 && continue
        H += pref*gg*( Qz[n]*Qz[mm] + 0.5*(Qp[n]*Qm[mm] + Qm[n]*Qp[mm]) )
    end

    # color-SINGLET projection (Gauss law on a ring ⇒ total color charge = 0):
    # the periodic Coulomb has no zero mode, so non-singlets are unpenalised —
    # lift them with Λ (Q_tot·Q_tot),  Q_tot^a = Σ_n Q_n^a.
    if Λ != 0
        Qzt = sum(Qz); Qpt = sum(Qp); Qmt = sum(Qm)
        H += Λ * ( Qzt*Qzt + 0.5*(Qpt*Qmt + Qmt*Qpt) )
    end
    return H, nq
end

# baryon number B = ¼ Σ σ^z ; return indices of basis states with B == 0
function b0_indices(nq::Int)
    idx = Int[]
    half = nq ÷ 2
    for s in 0:(2^nq-1)
        count_ones(s) == half && push!(idx, s+1)   # #up == #down  ⇒ Σσ^z = 0
    end
    return idx
end

# ---- quick test: vacuum + lowest excitations in B=0 -------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    for N in (4, 6)
        H, nq = build_su2_hamiltonian(N; a=1.0, g=5.0, m=0.2)
        idx = b0_indices(nq)
        Hb = Matrix(H[idx, idx])
        @assert norm(Hb - Hb') < 1e-9 "non-Hermitian"
        ev = sort(real.(eigvals(Hermitian(Hb))))
        gap = ev[2]-ev[1]
        println("N=$N (2N=$nq qubits, dim B0=$(length(idx))):")
        println("  E0=$(round(ev[1],digits=4))  gaps above vacuum: ",
                round.(ev[2:min(6,end)] .- ev[1], digits=4))
        println("  meson mass estimate 2m+3/8·ag² = ", 2*0.2 + (3/8)*25, "  (strong-coupling, OBC)")
    end
end
