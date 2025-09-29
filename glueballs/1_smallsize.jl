locops = (JLD2.load("glueballs/0_constants.jld2"))["single_stored_object"]
smallsize = (JLD2.load("glueballs/1_smallsize.jld2"))["single_stored_object"]
# smallsize = Dict()

groups = ["ZZ3", "SU3"]
lengths = [8, 9]
d = 3

# Compute disprel
if true
    for group in groups
    smallsize[group] = Dict()

    # Coupling loop
    λ_list = [1//10, 3//10, 5//10, 7//10, 9//10]

    for λ ∈ λ_list
    smallsize[group][λ] = Dict()

    # Small system size loop
    for l ∈ lengths
    smallsize[group][λ][l] = Dict()
    dct = smallsize[group][λ][l]

    # Local operators
    dct["h"] = λ * locops[group]["Esq3"] - (1 - λ) * (locops[group]["Up3"] + locops[group]["Up3"]')

    # Operators
    dct["C"] = kron_power(SparseMatrixCSC(locops[group]["C1"]), l)
    dct["T"] = operator_translation(SparseMatrixCSC, d, l)
    dct["P"] = operator_reflection(SparseMatrixCSC, d, l) # For this model is like that
    dct["H"] = summation_local(dct["h"], d, l; pbc = true)
    dct["R"] = dct["P"] * dct["C"]

    dct["disprel"] = Dict()

    # Dispersion relation
    drel = dispersion_relation(dct["H"], dct["T"], l, nlevels = 30)
    dct["disprel"]["states"] = drel
    dct["disprel"]["E"] = [energy(el) for el ∈ drel]
    dct["disprel"]["k"] = [momentum(el) for el ∈ drel]
    dct["disprel"]["C"] = [real(wavefunction(el)' * dct["C"] * wavefunction(el)) for el ∈ drel]

    # Get band properties
    gs = pop_groundstate!(drel)
    dct["groundstate"] = gs
    dct["gsenergy"] = energy(gs)
    firstband = get_firstband(drel)
    dct["firstband"] = firstband
    # We compute the average of the energies of the first band
    energiesband = []
    for el in firstband
        if el.koverpi == 0//1 || el.koverpi == 1//1 || el.koverpi == -1//1
            push!(energiesband, energy(el))
        else # We count two times the others
            push!(energiesband, energy(el))
            push!(energiesband, energy(el))
        end
    end
    dct["centroid"] = sum(energiesband) / length(energiesband)

    end # system size loop
    end # coupling loop
    end # group loop

    JLD2.save_object("glueballs/1_smallsize.jld2", smallsize)
end

function color_from_lambda(λ::Real)
    pE = [0.2,0.3,0.7]
    pB = [0.4,0.7,0.8]
    pI = λ * pE + (1 - λ) * pB
    return RGB(pI...)
end

# plot disprel
if true
    l = 8
    λ_list = [1//10, 3//10, 5//10, 7//10, 9//10]
    ncols = length(λ_list)
    for group in ["ZZ3", "SU3"]
        plt = Plots.plot(
            layout=(1,ncols), 
            legend = false,
            plot_padding = 0mm,
            margin = 0mm,
            framestyle = :box,
            guidefont = font(7),     # axis labels
            tickfont  = font(7),     # tick labels
            legendfont = font(7),    # legend
            titlefont = font(7)
            )
        for ii in 1:ncols
            λ = λ_list[ii]
            Egs = smallsize[group][λ][l]["gsenergy"]
            Ebar = smallsize[group][λ][l]["centroid"]
            ΔE = Ebar - Egs
            Elist = smallsize[group][λ][l]["disprel"]["E"]
            Klist = smallsize[group][λ][l]["disprel"]["k"]
            Clist = smallsize[group][λ][l]["disprel"]["C"]
            neg_C = findall(c -> c < 0, Clist)
            pos_C = findall(c -> c > 0, Clist)

            leb = L"\langle\mathcal{E}\rangle"
            ytickstoput = ([Egs, Egs + ΔE, Egs + 2ΔE, Egs + 3ΔE], [L"0", leb, L"2" * leb, L"3" * leb])
            leftmargin = 0mm
            if ii != 1
                ytickstoput = ([Egs, Egs + ΔE, Egs + 2ΔE, Egs + 3ΔE], ["", "", "", ""])
                leftmargin = -5mm
            end

            if group == "SU3"
                leb = L"\langle\mathcal{E}\rangle"
                ytickstoput = ([Egs, Egs + ΔE, Egs + 2ΔE], [L"0", leb, L"2" * leb])
                leftmargin = 0mm
                if ii != 1
                    ytickstoput = ([Egs, Egs + ΔE, Egs + 2ΔE], ["", "", ""])
                    leftmargin = -5mm
                end
            end

            titlelatex = L"\lambda = \dfrac{%$(numerator(λ))}{%$(denominator(λ))}"

            multipleformargin = 3
            if group == "SU3"
                multipleformargin = 2.5
            end

            Plots.scatter!(
                plt, 
                Klist[pos_C], 
                Elist[pos_C], 
                color = :blue, 
                subplot = ii, 
                markersize = 2, 
                xlims  = [0-0.3, π+0.3],
                xticks = ([0, π/2, π], [L"0", L"\dfrac{\pi}{2}", L"\pi"]),
                ylims  = [Egs - 0.05 * (3ΔE), Egs + multipleformargin * ΔE + 0.05 * (3ΔE)],
                yticks = ytickstoput,
                left_margin = leftmargin,
                title = titlelatex
                )
            Plots.scatter!(
                plt,
                Klist[neg_C], 
                Elist[neg_C], 
                color = :red, 
                marker = :diamond, 
                markersize = 2, 
                subplot = ii
                )
                print("Size of plt: ", size(Klist[pos_C]))
        end
        Plots.plot!(size=(220,180), grid=true)
        Plots.savefig("disprel$(group).pdf")
    end
end