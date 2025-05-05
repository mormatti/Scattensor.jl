using ITensors
using ITensorMPS
using Plots
using Colors
using Scattensor
using Revise
gr()

t = 25
Aelm = smatrix_expansion_real_space_tdvp(d, Hmpo, psi0big, W, χmax, t)
println("Done!")

function dofunc()
    # the size of one dimension of Aelm
    sampl = 100
    Kmatrix = zeros(ComplexF64, sampl, sampl)
    println("Computing the K matrix...")
    for i in 1:sampl
        for i′ in 1:sampl
            k = (i-1)/sampl * 2π - π
            k′ = (i′-1)/sampl * 2π - π
            Kmatrix[i,i′] = smatrix_element_momentum_space_tdvp(Aelm, π/2, -π/2, k, k′)
        end
    end
    # We eliminate the first 5 and last 5 rows and the first 5 and last 5 columns
    complex_colormap_plot(Kmatrix)
end
