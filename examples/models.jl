using Scattensor
using ITensors, ITensorMPS
using LinearAlgebra
using SparseArrays
using Optim
using Plots
using PlotlyJS
using KrylovKit
using Logging

# Hardcoding of local Hamiltonians of some 1D models.

function H0_spinchain(; Jx = 0, Jy = 0, Jz = 0, hx = 0, hy = 0, hz = 0, L0 = 2)
    os = OpSum()
    
    if L0 == 2
        os += -hx, "Sx", 1
        os += -hx, "Sx", 2
        os += -hy, "Sy", 1
        os += -hy, "Sy", 2
        os += -hz, "Sy", 1
        os += -hz, "Sy", 2
        os += -4*Jx, "Sx", 1, "Sx", 2
        os += -4*Jy, "Sy", 1, "Sy", 2
        os += -4*Jz, "Sz", 1, "Sz", 2
    elseif L0 == 3
        os += -2*hx, "Sx", 2
        os += -2*hy, "Sy", 2
        os += -2*hz, "Sy", 2
        os += -2*Jx, "Sx", 1, "Sx", 2
        os += -2*Jx, "Sx", 2, "Sx", 3
        os += -2*Jy, "Sy", 1, "Sy", 2
        os += -2*Jy, "Sy", 2, "Sy", 3
        os += -2*Jz, "Sz", 1, "Sz", 2
        os += -2*Jz, "Sz", 2, "Sz", 3
    else
        error("Local Hamiltonian is not defined for this number of sites.")
    end

    H0 = MPO(os, siteinds("S=1/2", L0))

    # We extract the site indices of the MPO
    sind = siteinds(H0)

    # We extract the indeces up and indeces down
    sind_trans = [collect(t) for t in zip(sind...)]
    sitesup = sind_trans[1]
    sitesdown = sind_trans[2]

    # We fuse the indices of the MPO to get a matrix
    # We transform the MPO into a matrix,
    # sequentially contract each MPO tensor into IT
    combinerup = combiner(sitesup)
    combinerdown = combiner(sitesdown)
    H0it = H0[1]
    for j in 2:L0
        H0it *= H0[j]
    end

    return SparseMatrixCSC(matrix(combinerup * H0it * combinerdown))
end

function H0_bosehubbard(; t = 0, U = 0, μ = 0, V = 0)
    os = OpSum()
    os += -t, "A", 1, "Adag", 2
    os += -t, "Adag", 1, "A", 2
    os += +U/4, "N", 1, "N", 1
    os += -U/4, "N", 1
    os += +U/4, "N", 2, "N", 2
    os += -U/4, "N", 2
    os += -μ/2, "N", 1, "N", 1
    os += -μ/2, "N", 2, "N", 2
    os += -V, "N", 1, "N", 2
    return MPO(os, siteinds("Boson", 2))
end

H0 = LocalOperator(H0_spinchain(Jz = 0.5, hx = 1, L0 = 3), [2, 2, 2], "h")

drel = disprel(H0, 11)

Plots.plot(drel)

gs = popgroundstate!(drel)

E0 = energy(gs)

bd = selectfirstband(drel)

energ, details = wannier(bd, H0, E0)

Plots.plot(energ)

print("Prova")