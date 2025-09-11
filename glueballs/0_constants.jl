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


locops = (JLD2.load("glueballs/locops.jld2"))["single_stored_object"]

Nlist = [2,3]
locops["ZZ"] = Dict()
locops["SU"] = Dict()
ZZ = locops["ZZ"]
SU = locops["SU"]

for N in Nlist
    ZZ[N] = Dict()
    SU[N] = Dict()

    # Local dimension
    ZZ[N]["N"] = N
    SU[N]["N"] = N
    ZZ[N]["localdim"] = N
    SU[N]["localdim"] = N

    # Identity matrix
    idm = operator_identity(SparseMatrixCSC, N)
    ZZ[N]["I1"] = idm
    ZZ[N]["I2"] = kron(idm, idm)
    ZZ[N]["I3"] = kron(idm, idm, idm)
    SU[N]["I1"] = idm
    SU[N]["I2"] = kron(idm, idm)
    SU[N]["I3"] = kron(idm, idm, idm)

    # Local pure reflection
    ZZ[N]["R2"] = operator_reflection(SparseMatrixCSC, N, 2)
    ZZ[N]["R3"] = operator_reflection(SparseMatrixCSC, N, 3)
    SU[N]["R2"] = operator_reflection(SparseMatrixCSC, N, 2)
    SU[N]["R3"] = operator_reflection(SparseMatrixCSC, N, 3)

    # The cyclic supradiagonal local operator
    tau = spzeros(N, N)
    for i in 1:N-1
        tau[i, i+1] = 1
    end
    tau[N, 1] = 1
    ZZ[N]["tau"] = tau
    SU[N]["tau"] = tau

    # The charge conjugation 1,2,3-local operator
    cc = spzeros(N, N)
    for i in 1:N
        cc[1 + mod(-(i-1), 0:N-1), i] = 1
    end
    ZZ[N]["C1"] = cc
    ZZ[N]["C2"] = kron(cc, cc)
    ZZ[N]["C3"] = kron(cc, cc, cc)
    SU[N]["C1"] = cc
    SU[N]["C2"] = kron(cc, cc)
    SU[N]["C3"] = kron(cc, cc, cc)

    kb(n,m) = ketbra(n,m,N)

    # Electric field (casimir) as a 2-local operator
    C2ZZ(N,n) = N^2 / π^2 * sin(π * n / N)^2
    C2SU(N,n) = (N + 1) // (2N) * (N - n) * n
    casimirZZN = spzeros(N^2, N^2)
    casimirSUN = spzeros(N^2, N^2)
    for n1 in 0:N-1
        for n2 in 0:N-1
            i1 = n1 + 1
            i2 = n2 + 1
            m = mod(n1 - n2, 0:N-1)
            E2ZZN = C2ZZ(N,n1) + C2ZZ(N,n2) + C2ZZ(N,m)
            E2SUN = C2SU(N,n1) + C2SU(N,n2) + C2SU(N,m)
            casimirZZN += E2ZZN * kron(kb(i1,i1), kb(i2,i2))
            casimirSUN += E2SUN * kron(kb(i1,i1), kb(i2,i2))
        end
    end
    ZZ[N]["Esq2"] = casimirZZN
    SU[N]["Esq2"] = casimirSUN

    # Electric field (casimir) as a 3-local operator
    ZZ[N]["Esq3"] = 1/2 * kron(ZZ[N]["Esq2"], ZZ[N]["I1"]) + 1/2 * kron(ZZ[N]["I1"], ZZ[N]["Esq2"])
    SU[N]["Esq3"] = 1/2 * kron(SU[N]["Esq2"], SU[N]["I1"]) + 1/2 * kron(SU[N]["I1"], SU[N]["Esq2"])

    # Plaquette term for ZN
    ZZ[N]["Up"] = ZZ[N]["tau"]
    ZZ[N]["Up3"] = kron(ZZ[N]["I1"], ZZ[N]["tau"], ZZ[N]["I1"])

    # Diagonal D terms for SU(N)
    DN = []
    if N == 2
        append!(DN, [Diagonal([1, -1/2])])
        append!(DN, [Diagonal([1, -1/2])])
    elseif N == 3
        append!(DN, [Diagonal([1, -1/√3, 1/3])])
        append!(DN, [Diagonal([1, 1/3, -1/√3])])
        append!(DN, [Diagonal([1, -1/√3, -1/√3])])
    elseif N == 4
        error("N = 4 is not implemented yet")
        append!(DN, [Diagonal([1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1])])
    elseif N == 5
        error("N = 5 is not implemented yet")
        append!(DN, [Diagonal([1, 1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1, 1])])
    elseif N == 6
        error("N = 6 is not implemented yet")
        append!(DN, [Diagonal([1, 1, 1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1, 1, 1])])
    elseif N == 7
        error("N = 7 is not implemented yet")
        append!(DN, [Diagonal([1, 1, 1, 1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1, 1, 1, 1])])
        append!(DN, [Diagonal([1, 1, 1, 1, 1, 1, 1])])
    end

    # Plaquette term for ZN
    Up = kron(DN[1], kb(N,1), DN[1])
    for n in 2:N
        Up += kron(DN[n], kb(n-1,n), DN[n])
    end
    SU[N]["Up"] = Up
    SU[N]["Up3"] = Up
end

JLD2.save_object("glueballs/locops.jld2", locops)