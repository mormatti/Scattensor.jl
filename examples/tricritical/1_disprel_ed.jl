
using SparseArrays
using Scattensor
using Revise
using Plots
using LinearAlgebra
using KrylovKit
using Printf

# --- Parameters ---
# --- Parameters ---
d = 3           # Local dimension (Spin-1)
L = 10          # System size (updated to 12)
L0 = 2          # Range of local Hamiltonian (nearest-neighbor)

ac = 0.91024
bc = 0.4156
gc = 0.0
theta = atan(2.224)
epsilon = 0
phi = 2.09

# Greek aliases for plotting
ϵ = epsilon
θ = theta
φ = phi

# Create plots directory if it doesn't exist
plots_dir = joinpath(@__DIR__, "plots")
mkpath(plots_dir)

println("Computing spectrum for L=$L")
println("Parameters: ϵ=$ϵ, θ=$θ, φ=$φ")

# --- Local Operators (Spin-1) ---
sz = spdiagm(0 => [1.0, 0.0, -1.0])
val = 1.0 / sqrt(2)
sx = spdiagm(1 => [val, val], -1 => [val, val])
id = spdiagm(0 => ones(d))

# Helper for tensor product
⊗(A, B) = kron(A, B)

# Squared operators
sx2 = sx * sx
sz2 = sz * sz

fname = "plot.png"
Hc = ac * (sx2 ⊗ id) + bc * (sz ⊗ id) + gc * (sz2 ⊗ id) - (sx ⊗ sx)
Hperp = -sin(theta) * (sx2 ⊗ id) + cos(theta) * (sz ⊗ id)
H0 = Hc + epsilon * Hperp

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

# --- Plotting ---
fname = "plot.png"
output_path = joinpath(plots_dir, fname)

title_str = @sprintf("L=%d, ϵ=%.1f, θ=%.2f, φ=%.2f", L, ϵ, θ, φ)
pl = plot_disprel(drel, -π, 2π; title=title_str)

savefig(pl, output_path)
println("Done. Plot saved to $output_path")
