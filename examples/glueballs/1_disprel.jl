using Scattensor
using ITensors, ITensorMPS
using LinearAlgebra
using SparseArrays
using Optim
using Plots
using PlotlyJS
using KrylovKit
using Logging
using Revise
plotlyjs()

# Constants
d = 3
L0 = 3
L = 9
λ = 0.5

# Local operators
c = [1 0 0; 0 0 1; 0 1 0]
E2 = Glueballs.electric_operator()
Up = Glueballs.plaquette_operator()
H0 = λ * E2 - (1 - λ) * (Up + Up')

# Operators
C = kron_power(SparseMatrixCSC(c), L)
T = operator_translation(SparseMatrixCSC, d, L)
P = operator_reflection(SparseMatrixCSC, d, L)
H = summation_local(H0, L0, d, L; pbc = true)
R = P * C

# Dispersion relation
drel = disprel(H, T, L, nlevels = 10)
plot_disprel(drel)