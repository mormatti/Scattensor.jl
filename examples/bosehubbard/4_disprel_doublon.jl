
using SparseArrays
using Scattensor
using Revise
using Plots
using LinearAlgebra
using KrylovKit
using Printf
using LaTeXStrings

# =====================================================================
#  Bose-Hubbard chain, d=3 (n_max=2), expecting band - continuum - band
# =====================================================================
#  H = -t Σ (b†_i b_{i+1} + h.c.) + (U/2) Σ n_i (n_i - 1) + ε Σ n_i
#  Large ε pushes vacuum to be GS and well-separates particle-number sectors.
#  For t=1, U=10, ε=20:
#    1-boson band:        [ε-2t, ε+2t]      = [18, 22]
#    2-boson continuum:   [2ε-4t, 2ε+4t]    = [36, 44]
#    Doublon band:        ~2ε + U ± O(t²/U) = 50 ± 0.2
#    3-boson continuum:   [3ε-6t, 3ε+6t]    = [54, 66]    → doublon below it
# =====================================================================

# --- Parameters ---
d = 3
L = 10
L0 = 2
t = 2.0
U = 10.0
ε = 20.0
nlev = 20

# --- Local operators ---
b   = SparseMatrixCSC(spdiagm(1 => [sqrt(i) for i in 1:(d-1)]))     # annihilation
bd  = SparseMatrixCSC(sparse(b'))                                    # creation
nop = SparseMatrixCSC(spdiagm(0 => Float64.(0:(d-1))))               # number
id3 = SparseMatrixCSC(spdiagm(0 => ones(d)))

⊗(A, B) = kron(A, B)

# --- Two-site local Hamiltonian ---
H0 = -t * (bd ⊗ b + b ⊗ bd)
H0 += 0.5 * (U/2) * (nop * (nop - id3)) ⊗ id3
H0 += 0.5 * (U/2) * id3 ⊗ (nop * (nop - id3))
H0 += 0.5 * ε * (nop ⊗ id3 + id3 ⊗ nop)

# --- Build H, T ---
println("Building H ...")
H = summation_local(H0, d, L; pbc = true)
println("dim(H) = $(size(H,1))")
println("Building T ...")
T = operator_translation(SparseMatrixCSC, d, L)

# --- Dispersion relation ---
println("Computing dispersion relation, nlevels = $nlev ...")
drel = dispersion_relation(H, T, L; nlevels = nlev)

# --- Spectrum table ---
println()
println("Spectrum (first $nlev levels per k):")
ks = sort(unique([bs.koverpi for bs in drel]))
println(rpad("k/π", 8), join([rpad(@sprintf("lev %d", i), 11) for i in 1:nlev]))
for kop in ks
    energies = sort([energy(bs) for bs in drel if bs.koverpi == kop])
    row = rpad(string(kop), 8)
    for i in 1:min(nlev, length(energies))
        row *= rpad(@sprintf("%.3f", energies[i]), 11)
    end
    println(row)
end

# --- Plot ---
plots_dir = joinpath(@__DIR__, "plots")
mkpath(plots_dir)
fname = @sprintf("disprel_bh_doublon_L%d_t%.1f_U%.1f_eps%.1f.png", L, t, U, ε)
output_path = joinpath(plots_dir, fname)

title_str = @sprintf("BH (d=3): L=%d, t=%.1f, U=%.1f, ε=%.1f", L, t, U, ε)
pl = plot_disprel(drel, -π, 2π;
    title = title_str, aspect_ratio = :auto, size = (900, 700))
savefig(pl, output_path)
println()
println("Plot saved to $output_path")
