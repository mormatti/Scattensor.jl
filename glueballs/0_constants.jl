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
using Measures
using LaTeXStrings

let
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

    locops = (JLD2.load("glueballs/0_constants.jld2"))["single_stored_object"]

    locops["ZZ3"] = Dict()
    locops["SU3"] = Dict()
    ZZ3 = locops["ZZ3"]
    SU3 = locops["SU3"]

    # Local dimension
    ZZ3["localdim"] = 3
    SU3["localdim"] = 3

    # Identity matrix
    idm = operator_identity(SparseMatrixCSC, 3)
    ZZ3["I1"] = idm
    ZZ3["I2"] = kron(idm, idm)
    ZZ3["I3"] = kron(idm, idm, idm)
    SU3["I1"] = idm
    SU3["I2"] = kron(idm, idm)
    SU3["I3"] = kron(idm, idm, idm)

    # Local pure reflection
    ZZ3["R2"] = operator_reflection(SparseMatrixCSC, 3, 2)
    ZZ3["R3"] = operator_reflection(SparseMatrixCSC, 3, 3)
    SU3["R2"] = operator_reflection(SparseMatrixCSC, 3, 2)
    SU3["R3"] = operator_reflection(SparseMatrixCSC, 3, 3)

    # The cyclic supradiagonal local operator
    tau = spzeros(3, 3)
    tau[1, 2] = 1
    tau[2, 3] = 1
    tau[3, 1] = 1
    ZZ3["tau"] = tau
    SU3["tau"] = tau

    # The charge conjugation 1,2,3-local operator
    cc = spzeros(3, 3)
    cc[1, 1] = 1
    cc[2, 3] = 1
    cc[3, 2] = 1
    ZZ3["C1"] = cc
    ZZ3["C2"] = kron(cc, cc)
    ZZ3["C3"] = kron(cc, cc, cc)
    SU3["C1"] = cc
    SU3["C2"] = kron(cc, cc)
    SU3["C3"] = kron(cc, cc, cc)

    kb(n,m) = ketbra(n,m,3)

    # Electric field (casimir) as a 2-local operator
    C2ZZ3(n) = 3^2 / π^2 * sin(π * n / 3)^2
    C2SU3(n) = 2 // 3 * (3 - n) * n
    casimirZZ3 = spzeros(3^2, 3^2)
    casimirSU3 = spzeros(3^2, 3^2)
    for n1 in 0:2
        for n2 in 0:2
            i1 = n1 + 1
            i2 = n2 + 1
            m = mod(n1 - n2, 0:2)
            E2ZZ3 = C2ZZ3(n1) + C2ZZ3(n2) + C2ZZ3(m)
            E2SU3 = C2SU3(n1) + C2SU3(n2) + C2SU3(m)
            casimirZZ3 += E2ZZ3 * kron(kb(i1,i1), kb(i2,i2))
            casimirSU3 += E2SU3 * kron(kb(i1,i1), kb(i2,i2))
        end
    end
    ZZ3["Esq2"] = casimirZZ3
    SU3["Esq2"] = casimirSU3

    # Electric field (casimir) as a 3-local operator
    ZZ3["Esq3"] = 1/2 * kron(ZZ3["Esq2"], ZZ3["I1"]) + 1/2 * kron(ZZ3["I1"], ZZ3["Esq2"])
    SU3["Esq3"] = 1/2 * kron(SU3["Esq2"], SU3["I1"]) + 1/2 * kron(SU3["I1"], SU3["Esq2"])

    # Plaquette term for ZN
    ZZ3["Up"] = ZZ3["tau"]
    ZZ3["Up3"] = kron(ZZ3["I1"], ZZ3["tau"], ZZ3["I1"])

    # Diagonal D terms for SU(N)
    Dg1 = Diagonal([1, 1/√3, 1/3])
    Dg2 = Diagonal([1, 1/3, 1/√3])
    Dg3 = Diagonal([1, 1/√3, 1/√3])

    # Plaquette term for ZN
    Up = kron(Dg1, kb(3,1), Dg1) + kron(Dg2, kb(1,2), Dg2) + kron(Dg3, kb(2,3), Dg3)

    SU3["Up"] = Up
    SU3["Up3"] = Up

    JLD2.save_object("glueballs/0_constants.jld2", locops)
end