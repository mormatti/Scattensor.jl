using ITensors
using ITensorMPS
using Plots
using Colors
using Scattensor
using Revise
gr()

N = 30
Aelm = smatrix_expansion_real_space(d, Hmpo, psi0big, W, N, χmax)
println("Done!")

function dofunc()
    Nmax = 30
    sampl = 100
    t = 35
    Kmatrix = zeros(ComplexF64, sampl, sampl)
    println("Computing the K matrix...")
    for n in [0,15,25,30]
        println("Step: $n")
        for i in 1:sampl
            for i′ in 1:sampl
                k = (i-1)/sampl * 2π - π
                k′ = (i′-1)/sampl * 2π - π
                Kmatrix[i,i′] = smatrix_element_momentum_space(Aelm, π/2, -π/2, k, k′, t; N = n)
            end
        end
        # We eliminate the first 5 and last 5 rows and the first 5 and last 5 columns
        complex_colormap_plot(Kmatrix)
        Plots.savefig("Glueb$n.png")
    end
end
