using Measures
using LaTeXStrings

locops = (JLD2.load("examples/glueballs/locops.jld2"))["single_stored_object"]
smallsize = (JLD2.load("examples/glueballs/smallsize.jld2"))["single_stored_object"]

groups = ["Z3"]

# Compute disprel
if true
    for group in groups
    
    d = locops[group]["localdim"]

    # Coupling loop
    λ_list = [0//1, 1//40, 1//30, 1//20, 1//10, 2//10, 3//10, 4//10, 5//10, 6//10, 7//10, 8//10, 9//10, 19//20, 29//30, 39//40, 1//1]
    for λ ∈ λ_list
    smallsize[group][λ] = Dict()

    # Small system size loop
    ls = [8]
    for l ∈ ls
    smallsize[group][λ][l] = Dict()
    dct = smallsize[group][λ][l]

    # Local operators
    dct["h"] = λ * locops[group]["Esq3"] - (1 - λ) * (locops[group]["Up3"] + locops[group]["Up3"]')

    # Operators
    dct["C"] = kron_power(SparseMatrixCSC(locops[group]["C1"]), l)
    dct["T"] = operator_translation(SparseMatrixCSC, d, l)
    dct["P"] = operator_reflection(SparseMatrixCSC, d, l)
    dct["H"] = summation_local(dct["h"], d, l; pbc = true)
    dct["R"] = dct["P"] * dct["C"]

    dct["disprel"] = Dict()

    # Dispersion relation
    drel = dispersion_relation(dct["H"], dct["T"], l, nlevels = 20)
    dct["disprel"]["E"] = [energy(el) for el ∈ drel]
    dct["disprel"]["k"] = [momentum(el) for el ∈ drel]
    dct["disprel"]["C"] = [real(wavefunction(el)' * dct["C"] * wavefunction(el)) for el ∈ drel]
    end # system size loop
    end # coupling loop

    end
end

JLD2.save_object("examples/glueballs/smallsize.jld2", smallsize)

function color_from_lambda(λ::Real)
    pE = [0.2,0.3,0.7]
    pB = [0.4,0.7,0.8]
    pI = λ * pE + (1 - λ) * pB
    return RGB(pI...)
end

# plot disprel
if true
    global sp = 0
    l = 8
    λ_list = [1//30, 1//10, 3//10, 5//10, 7//10, 9//10, 29//30]
    ncols = length(λ_list)
    p = Plots.plot(layout=(2,ncols), legend = false)
    for group in ["Z3", "SU(3)"]
        for ii in 1:ncols
            global sp = sp + 1
            print(sp)
            λ = λ_list[ii]
            Elist = smallsize[group][λ][l]["disprel"]["E"]
            Klist = smallsize[group][λ][l]["disprel"]["k"]
            Clist = smallsize[group][λ][l]["disprel"]["C"]
            neg_C = findall(c -> c < 0, Clist)
            pos_C = findall(c -> c > 0, Clist)
            Plots.scatter!(p, Klist[pos_C], Elist[pos_C], color = :blue, subplot = sp, markersize = 2, xlims = [0-0.3, π+0.3], xticks=([0, π/2, π], ["0", "π/2", "π"]), yticks=-100:100)
            Plots.scatter!(p, Klist[neg_C], Elist[neg_C], color = :red, marker = :diamond, markersize = 2, subplot = sp, xlims = [0-0.3, π+0.3])
        end
    end
    Plots.plot!(size=(700,400))
    Plots.savefig("AllLambdasPlots.pdf")
end