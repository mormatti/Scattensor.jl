
using SparseArrays
using Scattensor
using Revise
using Plots
using LinearAlgebra
using KrylovKit
using Printf

# =====================================================================
#  Effective dispersion relation from a trapped basis
# =====================================================================
#  Idea:
#   - H  = Σ_i h_i              (translation-invariant TFI)
#   - H' = Σ_i V(i) h_i         (trap: cos potential, PBC-friendly)
#   - Diagonalize H' → lowest N eigenstates |n⟩
#   - Build basis B = { T^p |n⟩ : n=1..N, p=0..L-1 }   (non-orthonormal)
#   - Compute overlap S, H_B, T_B in basis B
#   - Löwdin orthogonalize:  H̃ = S^{-1/2} H_B S^{-1/2}, T̃ = S^{-1/2} T_B S^{-1/2}
#   - Run dispersion_relation on (H̃, T̃, L)
# =====================================================================

# --- TFI parameters (same as 1_disprel_ed.jl) ---
d = 2
L = 10
L0 = 2
J = 1.0
h = 4.0

# --- Trap parameters ---
N = 5                  # number of seed eigenstates of H'
ic = L ÷ 2             # center of the trap (site index)
Vamp = 1.0             # amplitude of the cosine trap

# --- Local operators ---
σx = SparseMatrixCSC([0.0 1.0; 1.0 0.0])
σz = SparseMatrixCSC([1.0 0.0; 0.0 -1.0])
id2 = SparseMatrixCSC([1.0 0.0; 0.0 1.0])
⊗(A, B) = kron(A, B)

H0 = -J * (σz ⊗ σz) - 0.5 * h * (σx ⊗ id2) - 0.5 * h * (id2 ⊗ σx)

# --- Translation-invariant H (the one we will project later) ---
println("Constructing H (translation invariant)...")
H = summation_local(H0, d, L; pbc = true)
T = operator_translation(SparseMatrixCSC, d, L)

# --- Trap potential V(i) ---
#   V(i) = 1 - cos(2π (i - ic) / L)   ∈ [0, 2], min at i=ic, max at i=ic+L/2
V(i) = 1.0 - cos(2π * (i - ic) / L)

# --- Weighted Hamiltonian H' = Σ_j V(j) * T^j H0_ext T^{-j} -----------
println("Constructing H' (cosine trap)...")
Idext = operator_identity(SparseMatrixCSC, d^(L - L0))
H0ext = kron(H0, Idext)                          # H0 placed on sites (0,1)
Hprime = spzeros(ComplexF64, d^L, d^L)
Hcur = SparseMatrixCSC{ComplexF64}(H0ext)
for j in 0:(L-1)
    global Hcur, Hprime
    Hprime += V(j) * Hcur
    Hcur = T * Hcur * T'                          # shift to start at site j+1
end
Hprime = (Hprime + Hprime') / 2                   # symmetrize numerical noise

println("‖H' - H'†‖ = ", norm(Hprime - Hprime'))
println("‖[H', T]‖ = ", norm(Hprime * T - T * Hprime), "  (expected ≠ 0: trap breaks translations)")

# --- Diagonalize H', keep N lowest states ---
println("Diagonalizing H' for $N lowest levels...")
en_p, st_p, info = eigsolve(Hprime, size(Hprime, 1), N, :SR, ComplexF64;
                            ishermitian = true,
                            krylovdim = max(2N, N + 20))
en_p = real.(en_p[1:N])
st_p = st_p[1:N]
println("Lowest $N eigenvalues of H':")
for (i, e) in enumerate(en_p)
    @printf("  n=%2d   E_n' = %.8f\n", i, e)
end

# =====================================================================
#  Step 2 — build basis B = { T^p |n⟩ } and compute S, H_B, T_B
# =====================================================================
println()
println("Building basis B = {T^p |n⟩} (N=$N, L=$L → $(N*L) vectors)...")
dimH = size(H, 1)
Mdim = N * L
Vmat = zeros(ComplexF64, dimH, Mdim)
for n in 1:N
    Tpvn = Vector{ComplexF64}(st_p[n])
    for p in 0:(L-1)
        col = (n - 1) * L + p + 1
        Vmat[:, col] = Tpvn
        Tpvn = T * Tpvn
    end
end

println("Computing S = B†B, H_B = B†HB, T_B = B†TB ...")
S   = Vmat' * Vmat
H_B = Vmat' * (H * Vmat)
T_B = Vmat' * (T * Vmat)
S   = (S + S') / 2
H_B = (H_B + H_B') / 2

# =====================================================================
#  Step 3 — Löwdin orthogonalization with tolerance cutoff on S
# =====================================================================
println()
println("Löwdin orthogonalization (with cutoff on small S eigenvalues)...")
eS_vals, eS_vecs = eigen(Hermitian(S))
eS_vals = real.(eS_vals)
println("S eigenvalue range: [$(minimum(eS_vals)), $(maximum(eS_vals))]")
tol = 1e-8 * maximum(eS_vals)
keep = eS_vals .> tol
println("Keeping $(sum(keep)) / $Mdim modes (tol = $tol).")

eS_kept = eS_vals[keep]
US_kept = eS_vecs[:, keep]
W = US_kept * Diagonal(1 ./ sqrt.(eS_kept))         # M × Meff,  W† S W = I
Meff = size(W, 2)

Htilde = W' * H_B * W
Ttilde = W' * T_B * W
Htilde = (Htilde + Htilde') / 2

println("Effective subspace dimension Meff = $Meff")
println("‖H̃ - H̃†‖ = ", norm(Htilde - Htilde'))
println("‖T̃†T̃ - I‖ = ", norm(Ttilde' * Ttilde - I))
println("‖[H̃, T̃]‖ = ", norm(Htilde * Ttilde - Ttilde * Htilde))

# =====================================================================
#  Step 4 — dispersion relation from (H̃, T̃, L)
# =====================================================================
println()
println("Running dispersion_relation on the effective (H̃, T̃)...")
drel_eff = dispersion_relation(Matrix(Htilde), Matrix(Ttilde), L; nlevels = 10)

# Spectrum table
println()
println("Effective spectrum (first 10 levels per k):")
ks = sort(unique([bs.koverpi for bs in drel_eff]))
println(rpad("k/π", 10), join([rpad("level $i", 14) for i in 1:10]))
for kop in ks
    energies = sort([energy(bs) for bs in drel_eff if bs.koverpi == kop])
    row = rpad(string(kop), 10)
    for i in 1:min(10, length(energies))
        row *= rpad(@sprintf("%.6f", energies[i]), 14)
    end
    println(row)
end

# Plot
plots_dir = joinpath(@__DIR__, "plots")
mkpath(plots_dir)
fname = @sprintf("disprel_ising_eff_L%d_h%.2f_N%d.png", L, h, N)
output_path = joinpath(plots_dir, fname)
title_str = @sprintf("Effective spectrum (trap N=%d)  L=%d, h=%.2f", N, L, h)
pl = plot_disprel(drel_eff, -π, 2π;
                  title = title_str, aspect_ratio = :auto, size = (800, 600))
savefig(pl, output_path)
println("Plot saved to $output_path")
