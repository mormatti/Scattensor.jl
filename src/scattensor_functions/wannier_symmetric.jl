function wannier_symmetric(
    bandstates::Vector{<:BlochState}, # The states of the band,
    gsenergy::Real, # The groundstate energy of the system
    localhamop::Function, # The local hamiltonian of the system centered in the site of parity reflection
    parityop::Function, # The reflection operator of the system
    innerprod::Function # The inner product between Bloch state wavefunctions
    )

    info = Dict{String, Any}()

    println("Computing the maximally localized symmetric Wannier function.")

    # We check that all the momenta are expressed as fractions
    if !all(state -> state.koverpi isa Rational, bandstates)
        error("To select the momenta, the koverpi field must be a fraction.")
    end
    # We check that all momenta are positive
    if !all(state -> state.koverpi >= 0, bandstates)
        error("To get a symmetrized Wannier, only positive momenta are needed. Please select them.")
    end
    # We first check that ψ is really a single band, i.e. that two different states has different momentum.
    if !(length(unique(state.koverpi for state in bandstates)) == length(bandstates))
        error("ψ must be a single band.")
    end

    # We make lists of the momenta
    nonneg_ks = sort([state.koverpi for state in bandstates]) # All the non-negative momenta
    bulk_ks = sort(filter(x -> (x > 0 && x < 1), nonneg_ks)) # All the positive "bulk" momenta
    all_ks = sort(vcat(-bulk_ks, nonneg_ks)) # All the momenta (positive and negative)
    neg_ks = sort(-bulk_ks)

    # The length of the system is also the number of different momenta
    L = length(all_ks)

    # We create the dictionary k -> energy
    Engs = Dict(state.koverpi => energy(state) for state in bandstates)
    for k in neg_ks
        Engs[k] = Engs[-k]
    end

    # We compute the band centroid
    Ebar = sum(Engs[k] for k in all_ks) / L - gsenergy
    info["centroid"] = Ebar

    # We create the dictionary k -> wavefunction
    Ψ = Dict(state.koverpi => wavefunction(state) for state in bandstates)

    # We select the states with ψ = 0 and ψ = π, computing the parity of the band
    # If both exists, we check that they have the same parity
    # If only one exists, we take its parity
    # If none exists, we proceed with the symmetrization with parity = +1
    parity = 0
    parityguess(ψ1,ψ2) = Int64(round(real(innerprod(ψ1,parityop(ψ2)))))
    if haskey(Ψ, 0//1) && haskey(Ψ, 1//1)
        p0 = parityguess(Ψ[0//1],Ψ[0//1])
        pπ = parityguess(Ψ[1//1],Ψ[1//1])
        @assert p0 == pπ "The k = 0 and k = π states must have the same reflection symmetry."
        parity = p0
    elseif haskey(Ψ, 0//1)
        p0 = parityguess(Ψ[0//1],Ψ[0//1])
        parity = p0
    elseif haskey(Ψ, 1//1)
        pπ = parityguess(Ψ[1//1],Ψ[1//1])
        parity = pπ
    else
        parity = 1
    end
    print("Parity = $parity")
    println("")

    # We define the -k states from the +k states keeping the parity symmetry
    for k in neg_ks
        Ψ[k] = parity * parityop(Ψ[-k])
    end

    # We compute the effective band Hamiltonian Hc - E0/L as a matrix and a dictionary
    HcMinusGs = Dict((k1, k2) => innerprod(Ψ[k1], localhamop(Ψ[k2])) - (k1 == k2 ? gsenergy / L : 0) for k1 in all_ks, k2 in all_ks)

    # We compute the eigenvalues of the effective band Hamiltonian Hc - E0/L
    HcMatrix = [HcMinusGs[k1,k2] for k1 in all_ks, k2 in all_ks]
    eigHcMat = eigen(HcMatrix)
    println("Eigenvalues of local Hamiltonian: ")
    vals = eigHcMat.values
    info["Hcmatrix"] = HcMatrix
    info["Hceigvals"] = vals

    # We define the distance function
    χ(j) = sin(π * j / L)

    # We compute the effective band Hamiltonian Hj
    Hjeff = Dict((k1, k2, j) => exp(im * π * (k1 - k2) * j) * HcMinusGs[k1,k2] for k1 in all_ks, k2 in all_ks, j in -L:L)
    info["Hjeff"] = Hjeff

    # We define the spread operator normalized by L * Ebar
    # X2 = Dict((k1, k2) => sum(χ(j)^2 * Hjeff[k1,k2,j] / (L * Ebar) for j in 0:L-1) for k1 in all_ks, k2 in all_ks)

    # We define θ as the array of phases of non-negative ks minus the first
    param_ks = nonneg_ks[2:end]
    npars = length(param_ks)
    θ0 = rand(Float64, npars)

    function coeffs_from_thetas(θ)
        coeff = Dict{Rational{Int}, ComplexF64}()
        coeff[nonneg_ks[1]] = 1
        counter = 1
        for k in param_ks
            coeff[k] = exp(im * θ[counter])
            counter += 1
        end
        for k in neg_ks
            coeff[k] = coeff[-k]
        end
        return coeff
    end

    function funct(θ)
        coeff = coeffs_from_thetas(θ)
        num = 0
        den = 0
        for j in 0:L-1
            insideval = 0
            for k1 in all_ks
                for k2 in all_ks
                    insideval += coeff[k1]' * coeff[k2] * Hjeff[k1,k2,j]
                end
            end
            num += abs(χ(j)^2 * insideval)
            den += insideval
        end
        return real(num / den)
        # return Float64(real(sum(coeff[k1]' * coeff[k2] * X2[k1,k2] for k1 in all_ks for k2 in all_ks)))
    end

    # We find the optimal phases
    opt = optimize(funct, θ0, NelderMead())
    θopt = Optim.minimizer(opt)

    # We store the minimum
    info["minimum"] = Optim.minimum(opt)

    # We compute the phases and we save them
    phases = coeffs_from_thetas(θopt)
    info["phases"] = phases

    # We compute the density profile and we save it
    density = Dict(j => sum(real(phases[k1]' * phases[k2] * Hjeff[k1,k2,j]) for k1 in all_ks, k2 in all_ks) / (Ebar * L) for j in 0:Int(floor((L-1)/2)))
    info["density"] = density

    wopt = sum(phases[k] * Ψ[k] for k in all_ks)
    wopt = wopt / norm(wopt)

    print("Wannier function computed.")
    println("")

    return wopt, info
end

export wannier_symmetric