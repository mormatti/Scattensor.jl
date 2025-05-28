using ITensors
using ITensorMPS
using Plots
using Colors
using Scattensor
using Revise
gr()

# N = 30
t = 1
S = smatrix_real_space_tdvp(d, Hmpo, psi0big, W, t, nsteps, χmax)
println("S-matrix computed!")

function dofunc()
    Lk = 100
    matr = zeros(ComplexF64, Lk, Lk)
    println("Computing the matrix...")
    for i in 1:Lk
        for i′ in 1:Lk
            k = (i-1)/Lk * 2π - π
            k′ = (i′-1)/Lk * 2π - π
            matr[i,i′] = smatrix_element_momentum_space(S, k, -k, k′, -k′)
        end
    end
    # We eliminate the first 5 and last 5 rows and the first 5 and last 5 columns
    complex_colormap_plot(-matr)
    Plots.savefig("MatrixIsingTDVP,t=$t.png")
end
