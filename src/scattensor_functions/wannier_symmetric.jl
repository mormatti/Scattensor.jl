function wannier_symmetric(
    ψ::Vector{BlochState{Vector{<:Number}}}, # The states of the band,
    H0:: AbstractMatrix, # The local hamiltonian of the system
    L0::Integer, # The length of the support of the local hamiltonian
    H::AbstractMatrix, # The hamiltonian of the system
    L::Integer, # The length of the chain
    T::AbstractMatrix, # The translation operator of the system
    R::AbstractMatrix, # The reflection operator of the system
    d::Integer, # The local dimension of the system
    E0::Real, # The groundstate energy of the system
    )

    println("")
    print("Computing the maximally localized symmetric Wannier function.")
    println("")

    # TODO: WARNING: this function only with L odd and if k = 0 exists but not k = π
    # Fix it when it is possible

    # We first check that ψ is really a single band, i.e. that two different states has different momentum.
    length(unique(x.kfraction for x in ψ)) == length(ψ) || error("ψ must be a single band.")

    # We compute the sum of the energy of all the states
    H = H - I * E0
    H0ext = kron(H0, operator_identity(SparseMatrixCSC, d^(L-L0))) - E0 / L * I
    Etot = sum([(energy(ψ[α]) - E0) for α in eachindex(ψ)]) # We sum the energies of all the states

    Ψ = Dict(x.kfraction => wavefunction(x) for x in ψ) # We create a dictionary with the wavefunctions
    Engs = Dict(x.kfraction => energy(x) for x in ψ) # We create a dictionary with the wavefunctions
    all_keys = sort([ky for ky in keys(Ψ)])
    pos_keys = filter(x -> x > 0 && x < 1, keys(Ψ))

    speeds = []
    for i in 1:(length(all_keys)-1)
        mom(i) = all_keys[i]
        push!(speeds, 1/π * (Engs[mom(i+1)] - Engs[mom(i)]) / (mom(i+1) - mom(i)))
    end
    maxspeed = max(speeds...)

    print("Speeds: $speeds")
    println("")
    print("Maximum speed = $maxspeed")
    println("")

    # We select the states with ψ = 0 and ψ = π
    symmetrization_index = 0

    if haskey(Ψ, 0//1) && haskey(Ψ, 1//1)
        r0 = Int64(round(real(Ψ[0//1]' * R * Ψ[0//1])))
        rπ = Int64(round(real(Ψ[1//1]' * R * Ψ[1//1])))
        @assert r0 == rπ "The k=0 and k=π states must have the same reflection symmetry."
        symmetrization_index = r0
    elseif haskey(Ψ, 0//1)
        r0 = Int64(round(real(Ψ[0//1]' * R * Ψ[0//1])))
        symmetrization_index = r0
    elseif haskey(Ψ, 1//1)
        rπ = Int64(round(real(Ψ[1//1]' * R * Ψ[1//1])))
        symmetrization_index = rπ
    else
        # Proceeding with symmetrization
        symmetrization_index = 1
    end

    for kf in pos_keys
        rk = Ψ[-kf]' * R * Ψ[kf]
        Ψ[kf] = symmetrization_index / rk * Ψ[kf]
    end

    print("Symmetrization index: $symmetrization_index")
    println("")

    xc = (L-1)/2
    χ(x) = sin(π * (x - xc) / L) # We define the distance function.

    # We compute the matrix X2
    ψw = [Ψ[kf] for kf in all_keys]
    X2 = Hermitian([sum(χ(x)^2 * Ψ[k1]' * T^x * H0ext * (T')^x * Ψ[k2] for x in 0:L-1) for k1 in all_keys, k2 in all_keys])

    function funct(θ)
        w = exp.(im * [reverse(θ); 0; θ])
        ret = Float64((1/Etot) * real(w' * X2 * w))
        return ret # - Float64((1/trH)^2 * real(w' * X1 * w)^2)
    end

    opt = optimize(funct, rand(Float64, Int64((length(ψw)-1)/2)), SimulatedAnnealing())
    opt = optimize(funct, Optim.minimizer(opt), NelderMead())
    θopt = Optim.minimizer(opt)
    θopt = [reverse(θopt); 0; θopt]
    phasesopt = exp.(im * θopt)
    wopt = sum(phasesopt[α] * ψw[α] for α in eachindex(ψw))
    wopt = wopt / norm(wopt)

    print("Wannier function computed.")
    println("")

    return wopt
end

export wannier_symmetric