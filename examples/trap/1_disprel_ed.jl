
using SparseArrays
using Scattensor
using Revise
using Plots
using LinearAlgebra
using KrylovKit
using Printf

# --- Parameters ---
# Transverse-field Ising model:  H = -J Σ σz_i σz_{i+1} - h Σ σx_i
# Paramagnetic phase: h/J > 1 (transverse field dominates, Z2 symmetric).
d = 2           # Local dimension (spin-1/2)
L = 10          # System size
L0 = 2          # Range of local Hamiltonian (nearest-neighbor)

J = 1.0
h = 4.0         # h/J = 4 → deep inside the paramagnetic phase

# Create plots directory if it doesn't exist
plots_dir = joinpath(@__DIR__, "plots")
mkpath(plots_dir)

println("Computing spectrum for L=$L")
println("Parameters: J=$J, h=$h  (h/J = $(h/J), paramagnetic phase)")

# --- Local Operators (Pauli) ---
σx = SparseMatrixCSC([0.0 1.0; 1.0 0.0])
σz = SparseMatrixCSC([1.0 0.0; 0.0 -1.0])
id = SparseMatrixCSC([1.0 0.0; 0.0 1.0])

# Helper for tensor product
⊗(A, B) = kron(A, B)

# Two-site local Hamiltonian. The σx term is split symmetrically between the
# two sites so that summation_local with pbc reproduces -J ΣZZ - h ΣX exactly.
H0 = -J * (σz ⊗ σz) - 0.5 * h * (σx ⊗ id) - 0.5 * h * (id ⊗ σx)

# Summation to full Hamiltonian
println("Constructing Hamiltonian (SparseMatrixCSC)...")
H = summation_local(H0, d, L; pbc = true)
println("Hamiltonian size: $(size(H))")

# Translation Operator
println("Constructing Translation Operator for L=$L...")
T = operator_translation(SparseMatrixCSC, d, L)

# Compute Dispersion Relation / Spectrum
println("Computing Dispersion Relation (this may take a while)...")
drel = dispersion_relation(H, T, L; nlevels = 10)

# --- Print spectrum table: first 10 levels per momentum ---
println()
println("Spectrum (first 10 levels per k):")
ks = sort(unique([bs.koverpi for bs in drel]))
println(rpad("k/π", 10), join([rpad("level $i", 14) for i in 1:10]))
for kop in ks
    energies = sort([energy(bs) for bs in drel if bs.koverpi == kop])
    row = rpad(string(kop), 10)
    for i in 1:min(10, length(energies))
        row *= rpad(@sprintf("%.6f", energies[i]), 14)
    end
    println(row)
end
println()

# --- Plotting ---
fname = @sprintf("disprel_ising_para_L%d_h%.2f.png", L, h)
output_path = joinpath(plots_dir, fname)

title_str = @sprintf("TFI paramagnetic, L=%d, J=%.2f, h=%.2f", L, J, h)
pl = plot_disprel(drel, -π, 2π; title = title_str, aspect_ratio = :auto, size = (800, 600))

savefig(pl, output_path)
println("Done. Plot saved to $output_path")
