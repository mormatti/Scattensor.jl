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
using JLD2
plotlyjs()

locops = (JLD2.load("examples/glueballs/locops.jld2"))["single_stored_object"]

z2 = locops["Z2"]
z3 = locops["Z3"]
su2 = locops["SU(2)"]
su3 = locops["SU(3)"]

z2["localdim"] = 2
z3["localdim"] = 3
su2["localdim"] = 2
su3["localdim"] = 3

# Local operators
function ketbra(m::Integer, n::Integer, q::Integer)
    M = zeros(Int64, q, q)
    M[m, n] = 1
    return SparseMatrixCSC(M)
end

# Valid for both Z2 and SU(3)
function lambda_from_g(g::Real)
    return g^4 / (1 + g^4)
end

function g_from_lambda(λ::Real)
    if λ < 0 || λ >= 1
        error("Invalid λ value")
    end
    return (λ/(1 - λ))^(1/4)
end

# Local identity for Z2
z2["I1"] = operator_identity(SparseMatrixCSC, 2)
z2["I2"] = kron(z2["I1"], z2["I1"])
z2["I3"] = kron(z2["I1"], z2["I1"], z2["I1"])

# Local identity for Z3
z3["I1"] = operator_identity(SparseMatrixCSC, 3)
z3["I2"] = kron(z3["I1"], z3["I1"])
z3["I3"] = kron(z3["I1"], z3["I1"], z3["I1"])

# Local identity for SU(2)
su2["I1"] = operator_identity(SparseMatrixCSC, 2)
su2["I2"] = kron(su2["I1"], su2["I1"])
su2["I3"] = kron(su2["I1"], su2["I1"], su2["I1"])

# Local identity for SU(3)
su3["I1"] = operator_identity(SparseMatrixCSC, 3)
su3["I2"] = kron(su3["I1"], su3["I1"])
su3["I3"] = kron(su3["I1"], su3["I1"], su3["I1"])

# Local pure reflection for Z2
z2["R2"] = operator_reflection(SparseMatrixCSC, 2, 2)
z2["R3"] = operator_reflection(SparseMatrixCSC, 2, 3)

# Local pure reflection for Z3
z3["R2"] = operator_reflection(SparseMatrixCSC, 3, 2)
z3["R3"] = operator_reflection(SparseMatrixCSC, 3, 3)

# Local pure reflection for SU(2)
su2["R2"] = operator_reflection(SparseMatrixCSC, 2, 2)
su2["R3"] = operator_reflection(SparseMatrixCSC, 2, 3)

# Local pure reflection for SU(3)
su3["R2"] = operator_reflection(SparseMatrixCSC, 3, 2)
su3["R3"] = operator_reflection(SparseMatrixCSC, 3, 3)

# Charge conjugation for Z2
z2["C1"] = SparseMatrixCSC([1 0; 0 1])
z2["C2"] = kron(z2["C1"],z2["C1"])
z2["C3"] = kron(z2["C1"],z2["C1"],z2["C1"])

# Charge conjugation for Z3
z3["C1"] = SparseMatrixCSC([1 0 0; 0 0 1; 0 1 0])
z3["C2"] = kron(z3["C1"],z3["C1"])
z3["C3"] = kron(z3["C1"],z3["C1"],z3["C1"])

# Charge conjugation for SU(2)
su2["C1"] = SparseMatrixCSC([1 0; 0 1])
su2["C2"] = kron(su2["C1"],su2["C1"])
su2["C3"] = kron(su2["C1"],su2["C1"],su2["C1"])

# Charge conjugation for SU(3)
su3["C1"] = SparseMatrixCSC([1 0 0; 0 0 1; 0 1 0])
su3["C2"] = kron(su3["C1"],su3["C1"])
su3["C3"] = kron(su3["C1"],su3["C1"],su3["C1"])

# X reflection operator for Z2
z2["X2"] = z2["R2"]
z2["X3"] = z2["R3"]

# X reflection operator for Z3
z3["X2"] = z3["R2"] * z3["C2"]
z3["X3"] = z3["R3"] * z3["C3"]

# X reflection operator for SU(2)
su2["X2"] = su2["R2"]
su2["X3"] = su2["R3"]

# X reflection operator for SU(3)
su3["X2"] = su3["R2"] * su3["C2"]
su3["X3"] = su3["R3"] * su3["C3"]

# Y reflection operator for Z2
z2["Y2"] = z2["I2"]
z2["Y3"] = z2["I3"]

# Y reflection operator for Z3
z3["Y2"] = z3["C2"]
z3["Y3"] = z3["C3"]

# Y reflection operator for SU(2)
su2["Y2"] = su2["I2"]
su2["Y3"] = su2["I3"]

# Y reflection operator for SU(3)
su3["Y2"] = su3["C2"]
su3["Y3"] = su3["C3"]

# Tau operator for Z2
kb(n,m) = ketbra(n,m,2)
z2["τ"] = kb(1,2) + kb(2,1)

# Tau operator for Z3
kb(n,m) = ketbra(n,m,3)
z3["τ"] = kb(3,1) + kb(1,2) + kb(2,3)

# Tau operator for SU(2)
kb(n,m) = ketbra(n,m,2)
su2["τ"] = kb(1,2) + kb(2,1)

# Tau operator for SU(3)
kb(n,m) = ketbra(n,m,3)
su3["τ"] = kb(3,1) + kb(1,2) + kb(2,3)

# Electric field squared for Z2
z2["Esq2"] = 4/(π^2) * Diagonal([0,2,2,2])

# Electric field squared for Z3
z3["Esq2"] = 27/(2 * π^2) * Diagonal([0,2,2,2,2,3,2,3,2])

# Electric field squared for Z2
su2["Esq2"] = 3/4 * Diagonal([0,2,2,2])

# Electric field squared for SU(3)
su3["Esq2"] = 4/3 * Diagonal([0,2,2,2,2,3,2,3,2])

# Electric field squared as 3-local operator
z2["Esq3"] = 1/2 * kron(z2["Esq2"], z2["I1"]) + 1/2 * kron(z2["I1"], z2["Esq2"])
z3["Esq3"] = 1/2 * kron(z3["Esq2"], z3["I1"]) + 1/2 * kron(z3["I1"], z3["Esq2"])
su2["Esq3"] = 1/2 * kron(su2["Esq2"], su2["I1"]) + 1/2 * kron(su2["I1"], su2["Esq2"])
su3["Esq3"] = 1/2 * kron(su3["Esq2"], su3["I1"]) + 1/2 * kron(su3["I1"], su3["Esq2"])

# Plaquette term for Z2
z2["Up"] = z2["τ"]
z2["Up3"] = kron(z2["I1"], z2["τ"], z2["I1"])

# Plaquette term for Z3
z3["Up"] = z3["τ"]
z3["Up3"] = kron(z3["I1"], z3["τ"], z3["I1"])

# Plaquette term for SU(2)
su2["Up"] = SparseMatrixCSC(kron([1/2 0; 0 -1/4], [0 1; 1 0], [1/2 0; 0 -1/4]))
su2["Up3"] = su2["Up"]

# Plaquette term for SU(3)
kb(n,m) = ketbra(n,m,9)
su3["V"] = kb(1,1) + kb(3,2) + kb(2,3) + kb(5,4) + kb(4,5) + kb(6,6) + kb(9,7) + kb(8,8) + kb(7,9)
# su3["K"] = Diagonal([1, 1, 1, 3^(-1/4), 3^(-1/2), 3^(-1/4), 3^(-1/2), 3^(-1/4), 3^(-1/4)])
su3["K"] = Diagonal([1, 1, 1, -3^(-1/4), 3^(-1/2), -3^(-1/4), 3^(-1/2), -3^(-1/4), -3^(-1/4)])
su3["Kp"] = su3["V"] * su3["K"] * su3["V"]
su3["XY"] = su3["X2"] * su3["Y2"]
su3["XY3"] = su3["X3"] * su3["Y3"]
su3["KTL"] = su3["K"]
su3["KBL"] = su3["Y2"] * su3["Kp"] * su3["Y2"]
su3["KBR"] = su3["XY"] * su3["K"] * su3["XY"]
su3["KTR"] = su3["X2"] * su3["Kp"] * su3["X2"]
su3["Up"] = kron(su3["I1"], su3["τ"], su3["I1"]) * kron(su3["KTL"], su3["I1"]) * kron(su3["KBL"], su3["I1"]) * kron(su3["I1"], su3["KBR"]) * kron(su3["I1"], su3["KTR"])
su3["Up3"] = su3["Up"]

JLD2.save_object("examples/glueballs/locops.jld2", locops)