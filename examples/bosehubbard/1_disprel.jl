using SparseArrays
using Scattensor
using Revise
using Plots
using LinearAlgebra
using ITensorMPS, ITensors
using PlotlyJS
plotly()

d = 3
L = 9
L0 = 3

t = 0.1  # 0.15
U = 0   # 0
mu = 5/6   # 1
V = (-9-4*√2)/54  # 0.5

# Local operators
id = spdiagm(0 => ones(d))
b = spdiagm(-1 => [√i for i in 1:(d-1)])
n = spdiagm(0 => [i-1 for i in 1:d])

⊗(A, B) = kron(A, B)
H0 = -t/2 * (b ⊗ b' ⊗ id + id ⊗ b ⊗ b')
H0 += -t/2 * (b' ⊗ b ⊗ id + id ⊗ b' ⊗ b)
H0 += U/2 * (id ⊗ (n * (n - id)) ⊗ id)
H0 += -mu * (id ⊗ n ⊗ id)
H0 += V/2 * (id ⊗ n ⊗ n + n ⊗ n ⊗ id)

# Operators
T = operator_translation(SparseMatrixCSC, d, L)
R = operator_reflection(SparseMatrixCSC, d, L)
H = summation_local(H0, L0, d, L; pbc = true)

# Dispersion relation
drel = disprel(H, T, L, nlevels = 10)
plot_disprel(drel)
