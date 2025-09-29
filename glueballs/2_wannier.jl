let
    constants = (JLD2.load("glueballs/0_constants.jld2"))["single_stored_object"]
    smallsize = (JLD2.load("glueballs/1_smallsize.jld2"))["single_stored_object"]
    # wannier = (JLD2.load("glueballs/2_wannier.jld2"))["single_stored_object"]
    wannier = Dict()

    groups = ["ZZ3", "SU3"]
    lambdas = [1//10, 3//10, 5//10, 7//10, 9//10]
    lengths = [9]
    L0 = 3

    if true
        d = 3

        for group in groups
            wannier[group] = Dict()

            for λ ∈ lambdas
                wannier[group][λ] = Dict()

                for Li in lengths
                    wannier[group][λ][Li] = Dict()
                
                    dct0 = constants[group]
                    dct1 = smallsize[group][λ][Li]
                    dct2 = wannier[group][λ][Li]

                    cutoff = 1e-10
                    dct2["cutoff"] = cutoff
                    drel = dct1["disprel"]["states"] # Be careful, define this in the previous file
                    bpsi0 = pop_groundstate!(drel)
                    psi0 = mps_from_vector(wavefunction(bpsi0), d; cutoff = cutoff)
                    E0 = dct1["gsenergy"]
                    firstband = get_firstband(drel)
                    Hi = dct1["H"]
                    Ti = dct1["T"]
                    Ri = dct1["R"]
                    H0 = dct1["h"]
                    print("Firstband type = ", typeof(firstband), "\n")
                    psiw = wannier_symmetric(firstband, H0, L0, Hi, Li, Ti, Ri, d, E0)
                    density = local_exp_value(H0, psiw, L0, d, Li, addconst = -E0/Li)
                    centroid = dct1["centroid"]
                    densitynormalized = density / centroid
                end
            end
        end
    end

    if true
        λ = 3//10
        Li = 8
        group = "ZZ3"
        H0 = smallsize[group][λ][Li]["h"]
        psiw = wannier[group][λ][Li]["psiw"]
        E0 = smallsize[group][λ][Li]["gsenergy"]
        L0 = 3
        d = 3
        Plots.plot(log.(abs.(local_exp_value(H0, psiw, L0, d, Li, addconst = -E0/Li))))
    end

end