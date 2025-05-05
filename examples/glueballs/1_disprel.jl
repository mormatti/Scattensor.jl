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
Li = 9
λ = 0.5

# Local operators
c = [1 0 0; 0 0 1; 0 1 0]
E2 = Glueballs.electric_operator()
Up = Glueballs.plaquette_operator()
H0 = λ * E2 - (1 - λ) * (Up + Up')

# Operators
Ci = kron_power(SparseMatrixCSC(c), Li)
Ti = operator_translation(SparseMatrixCSC, d, Li)
Pi = operator_reflection(SparseMatrixCSC, d, Li)
Hi = summation_local(H0, L0, d, Li; pbc = true)
Ri = Pi * Ci

# Dispersion relation
drel = disprel(Hi, Ti, Li, nlevels = 10)
p = plot_disprel(drel, color = hex_to_rgb("#CD6C2E"))