using SparseArrays
using Scattensor
using Revise
using Plots
using LinearAlgebra
using ITensorMPS, ITensors

# Constants
d = 2
L0 = 3
L = 11
λ = 0.2

# Local operators
σx = [0 1; 1 0]
σz = [1 0; 0 -1]
Id = [1 0; 0 1]
H0 = SparseMatrixCSC(-λ * (1/2 * kron(Id, σz, σz) + 1/2 * kron(σz, σz, Id)) - (1 - λ) * (kron(Id, σx, Id)))

# Operators
H = summation_local(H0, L0, d, L; pbc = true)
T = operator_translation(SparseMatrixCSC, d, L)
R = operator_reflection(SparseMatrixCSC, d, L)

drel = disprel(H, T, L, nlevels = 10)

plot_disprel(drel)