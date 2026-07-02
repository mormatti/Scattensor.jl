
using SparseArrays
using Scattensor
using Revise
using Plots
using LinearAlgebra
using ITensors, ITensorMPS

# --- Parameters ---
d = 2           # Local dimension (Spin-1/2)
L = 10          # System size (Tensor Network method)
J = 0.2         # Interaction strength
g = 1.0         # Transverse field (paramagnetic-like phase)

println("Ising Model (TN Method)")
println("L = $L, J = $J, g = $g")

# --- Local Operators ---
σz = spdiagm(0 => [1.0, -1.0])
σx = spdiagm(1 => [1.0], -1 => [1.0])
id = spdiagm(0 => ones(d))

# Helper for tensor product
⊗(A, B) = kron(A, B)

# --- Local Hamiltonian Construction ---
# H = -J ∑ σz_i σz_{i+1} - g ∑ σx_i
term_interaction = -J * (σz ⊗ σz)
term_field = -g * (σx ⊗ id)

H0_matrix = term_interaction + term_field

# --- Convert to MPO and Sum ---
println("Converting local Hamiltonian to MPO...")
H0_mpo = mpo_from_matrix(H0_matrix, d)

println("Summing local MPO over the chain...")
H = summation_local(H0_mpo, L; pbc = true)

# --- Dispersion Relation (TN/DMRG) ---
println("Computing Dispersion Relation (nlevels=5, maxdim=200, tol=1e-6)...")
drel = dispersion_relation(H; nlevels = 5, maxdim = 200, tol = 1e-6)

# --- Plotting ---
println("Plotting...")
# Plot from -π to 2π
pl = plot_disprel(drel, -π, 2π)
savefig(pl, joinpath(@__DIR__, "disprel_ising_tn_L10_J02_g10.png"))
println("Done. Plot saved to $(joinpath(@__DIR__, "disprel_ising_tn_L10_J02_g10.png"))")
