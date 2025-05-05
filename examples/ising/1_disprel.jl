using SparseArrays
using Scattensor
using Revise
using Plots
using LinearAlgebra
using ITensorMPS, ITensors

# Constants
d = 2
L0 = 3
Li = 11
λ = 0.5

# Local operators
id = SparseMatrixCSC([1 0; 0 1])
σx = SparseMatrixCSC([0 1; 1 0])
σz = SparseMatrixCSC([1 0; 0 -1])
hunit = 4 / π^2
ME = hunit * spdiagm(0 => [2, 2, 2, 0])

⊗(A, B) = kron(A, B)
HE = 1/2 * (ME ⊗ id' + id ⊗ ME)
HB = id ⊗ σx ⊗ id + id ⊗ σx' ⊗ id
H0 = SparseMatrixCSC(λ * HE - (1 - λ) * HB)

# Operators
Hi = summation_local(H0, L0, d, Li; pbc = true)
Ti = operator_translation(SparseMatrixCSC, d, Li)
Ri = operator_reflection(SparseMatrixCSC, d, Li)

drel = disprel(Hi, Ti, Li, nlevels = 10)

plot_disprel(drel)