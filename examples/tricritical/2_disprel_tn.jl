
using SparseArrays
using Scattensor
using Revise
using Plots
using LinearAlgebra
using ITensors, ITensorMPS

# --- Parameters ---
d = 3           # Local dimension (Spin-1)
L = 30          # System size (Tensor Network method)
L0 = 2          # Range of local Hamiltonian (nearest-neighbor)

println("Blume-Capel Model (TN Method)")
println("L = $L")

# Tricritical point parameters
J_tric = 1.0
D_tric = 0.91024
Γ_tric = 0.41563
h_tric = 0.0

# Perturbation vector Δp = (ΔJ, ΔD, ΔΓ, Δh)
# Using factor 5 perturbation from previous step
# Perturbation vector Δp = (ΔJ, ΔD, ΔΓ, Δh)
strength = 5
ΔJ, ΔD, ΔΓ, Δh = 0.0, 0.0, -0.015*strength, 0.015*strength

J = J_tric + ΔJ
D = D_tric + ΔD
Γ = Γ_tric + ΔΓ
h = h_tric + Δh

println("Hamiltonian Parameters (with perturbation):")
println("J = $J, D = $D, Γ = $Γ, h = $h")

# --- Local Operators (Spin-1) ---
sz = spdiagm(0 => [1.0, 0.0, -1.0])
val = 1.0 / sqrt(2)
sx = spdiagm(1 => [val, val], -1 => [val, val])
id = spdiagm(0 => ones(d))

# Helper for tensor product
⊗(A, B) = kron(A, B)

# --- Local Hamiltonian Matrix Construction ---
# H(J,D,Γ,h) = -J ∑ S_i^z S_{i+1}^z + D ∑ (S_i^z)^2 - Γ ∑ S_i^x - h ∑ S_i^z

# Term 1: Interaction -J sz ⊗ sz
term_interaction = -J * (sz ⊗ sz)

# Term 2: Single-ion D (sz)^2
sz2 = sz * sz
term_D = D * (sz2 ⊗ id)

# Term 3: Transverse field -Γ sx
term_Γ = -Γ * (sx ⊗ id)

# Term 4: Longitudinal field -h sz
term_h = -h * (sz ⊗ id)

H0_matrix = term_interaction + term_D + term_Γ + term_h

# --- Convert to MPO and Sum ---
println("Converting local Hamiltonian to MPO...")
H0_mpo = mpo_from_matrix(H0_matrix, d)

println("Summing local MPO over the chain...")
H = summation_local(H0_mpo, L; pbc = true)

# --- Dispersion Relation (TN/DMRG) ---
println("Computing Dispersion Relation (nlevels=3, maxdim=50)...")
# Note: TN dispersion relation is more expensive, using nlevels=3
drel = dispersion_relation(H; nlevels = 3, maxdim = 50)

# --- Plotting ---
println("Plotting...")
# Plot from -π to π
pl = plot_disprel(drel, -π, π)
savefig(pl, joinpath(@__DIR__, "disprel_bc_tn_L30.png"))
println("Done. Plot saved to $(joinpath(@__DIR__, "disprel_bc_tn_L30.png"))")
