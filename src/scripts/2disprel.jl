include("models.jl")

""" Computes the dispersion relation of a system with a local Hamiltonian `H0`.

    # Assumptions

    - 1D: the system is a 1D chain of sites with unit spacing a = 1.

    - Uniformity: the system is assumed to be uniform, i.e. the local dimension d is 
      the same for all sites.

    - Closure: the system is closed and the Hamiltonian is Hermitian.

    - Translational invariance: the Hamiltonian is assumed to be invariant under 
      translations.

    - Locality: the Hamiltonian is local, i.e. it can be decompose in a sum of terms,
      which are the same but shifted by a certain number of sites, whose length
      do not scale with the system size.

    # Inputs (`args`):

    The input are all keyboards arguments.

    You have to pass:

    - The local Hamiltonian, which can be passed as:

        - `H0` in an MPO form (possibly sparse if you want to use the Krylov method).
        
        or:

        - `H0` in a matrix form (possibly sparse if you want to use the Krylov method).
        - `L0` the number of sites of the local Hamiltonian.
    
    You should pass (they are important arguments):

    - `Lrange`: the number of sites of the chain for which the dispersion relation 
       is calculated. It can be also a list if you want to calculate the dispersion
       relation for multiple lengths.
       Default is L0.

    - `nlevels` the number of levels to consider for each momentum.
       Default is 1.

    - `method` a symbol indicating the method. Default is `:krylov`, but can be 
       also `:exact`. Soon it will be available also `:tns`, using pure tensor networks.
    
    You can pass (facultative arguments):

    - `params` a dictionary containing the parameters of the model. It can be used
       to identify the model.

    # Output:

    A dictionary (or list of dictionaries) containing the following keys:

    - `energies::Dict`, a dictionary {momentum => energies::Real}.

    - `states::Dict`, a dictionary {momentum => states}.

    - `H0::Int`, see above.

    - `L0::int`, see above.

    - `d::Int` the local dimension of the system.

    - `method::Symbol` the method used to compute the dispersion relation.

    - `params::Dict` see above.
    """
function disprel(args...)
    # Managing the inputs

    # Keyword argument "Lrange"
    if haskey(kwargs, :Lrange)
        Lrange = kwargs[:Lrange]
    else
        Lrange = [L0]
    end
    
    # Keyword argument "method"
    method = :krylov
    if haskey(kwargs, :method)
        method = kwargs[:method]
        if method == :tns
            error("Pure Tensor Network method still not available.")
        elseif method == :exact
            @info "Proceeding with exact diagonalization."
        elseif method == :krylov
            @info "Proceeding with krylov decomposition."
        else
            error("Method not valid. Try with :exact or :krylov.")
        end
    else
        println("Method not specified, proceeding with krylov.")
    end
    otherdata[:method] = method

    # Keyword argument Lrange
    nlevels = 1
    if haskey(kwargs, :nlevels)
        nlevels = kwargs[:nlevels]
    end

    # First case: the input local Hamiltonian is a matrix
    if H0 isa AbstractMatrix
        # The matrix H0 in this case is the input local Hamiltonian
        H0mat = H0

        #assert that the matrix is squared
        N, N2 = size(H0mat)
        @assert N == N2

        # We extract the local dimension
        if haskey(kwargs, :L0)
            L0 = kwargs[:L0]
            d = Int(round(N^(1/L0)))
        else
            error("Local Hamiltonian system size must be provided.")
        end

    # Second case: the input local Hamiltonian is a MPO
    elseif H0 isa MPO
        # We extract the site indices of the MPO
        sind = siteinds(H0)

        # We extract the indeces up and indeces down
        sind_trans = [collect(t) for t in zip(sind...)]
        sitesup = sind_trans[1]
        sitesdown = sind_trans[2]

        # We extract the size of the local Hamiltonian
        L0 = length(sind)
        println("Support of Local Hamiltonian: L₀ = $L0. ")

        # We extract the local dimension
        d = dim(sind[1][1])
        println("Local dimension d = $d. ")

        # We check that the MPO is uniform in dimension
        for i in 1:L0
            for j in 1:2
                if d != dim(sind[i][j])
                    error("Wrong site dimensions of the MPO, it must be uniform.")
                end
            end
        end
        print("MPO is uniform.")

        # We fuse the indices of the MPO to get a matrix
        combinerup = combiner(sitesup)
        combinerdown = combiner(sitesdown)

        # We transform the MPO into a matrix,
        # sequentially contract each MPO tensor into IT
        H0it = H0[1]
        for j in 2:L0
            H0it *= H0[j]
        end

        # We convert the ITensor into a matrix 
        H0mat = matrix(combinerup * H0it * combinerdown)
    else
        error("Non valid input type.")
    end

    # We cut all the term inside L0 wich are smaller than L0
    Lrange = [L for L in Lrange if L >= L0]

    for L in Lrange

        # Def: the number of sites to add
        sc = L - L0

        # Def: the identity matrix for the extended chain (length L)
        Idm = spzeros(Float64, d^L, d^L)

        # We construct Idm
        for i in 1:d^L
            Idm[i, i] = 1.0
        end

        # Def: the identity matrix for the complementary chain (length L - L0)
        Idmext = spzeros(Float64, d^sc, d^sc)

        # We construct Idmext
        for i in 1:d^sc
            Idmext[i, i] = 1.0
        end
        
        # Def: the local Hamiltonian (at site 0) written for the extended chain
        H0matext = kron(H0mat, Idmext)

        # Def: the translation operator, which is a Sparse Matrix by construction
        T = generate_translation_operator_matrix(d, L)

        # Def: the Hamiltonian for the extended chain
        H = deepcopy(H0matext)

        # We construct H
        for _ in 1:L-1
            H0matext = T * H0matext * T'
            H += H0matext
        end

        # Def: the highest and lowest energy of the Hamiltonian
        Emax = 0.0
        Emin = 0.0

        # We extract the highest and lowest energy of the Hamiltonian
        if method == :krylov
            @info "Computing the groundstate energy..."
            val, _, _ = eigsolve(H, size(H, 1), 1, :SR, eltype(H); ishermitian = true)
            Emin = real(val[1])
            display(val)

            @info "Computing the highest energy..."
            val, _, _ = eigsolve(H, size(H, 1), 1, :LR, eltype(H); ishermitian = true)

            Emax = real(val[1])

        elseif method == :exact
            # We convert the Hamiltonian into a Hermitian matrix, since Hermitian
            # is a dense type
            H = Hermitian(H)

            @info "Computing the highest energy..."
            energies, _ = eigen(H)

            Emax = maximum(energies)
            Emin = minimum(energies)
        end

        # Def: the energy range.
        ΔE = Emax - Emin

        # Def: the energy range divisor (see later)
        div_fraction = 10

        # Def: the energy shift
        λ = Emax + ΔE / div_fraction

        # We project the Hamiltonian to the eigenspace of every momentum k and we
        # compute the dispersion relation.

        # Def: the vectors of energies
        energies = Dict{Rational, Vector{Real}}()

        # Def: the vectors of states
        states = Dict{Rational, Vector{Vector{Complex}}}()

        # Def: the projector to the momentum k (generating function)
        function P(k)
            proj = Idm
            for _ in 1:L-1
                proj = exp(-im * k) * T * proj + Idm
            end
            return proj
        end

        @info "Projecting and diagonalizing the Hamiltonian for each momentum k..."
        print("Momentum: ")
        # We iterate over all the possible momenta
        for kl in 0:(L-1)

            # We print the momentum
            print("$(kl)π ")

            # Def: f = k / π (which is assumed to be a fraction)
            f = (2 * kl) // L
            if f > 1
                f -= 2
            end
            k = π * f

            # Def: the Hamiltonian to project to the momentum k
            Hk = deepcopy(H)    

            # Only because the Hamiltonian is translational invariant,
            # we can apply the projector only to the left
            Hk = Hk - λ * Idm
            Hk = (P(k) * Hk) / L
            Hk =  Hk + λ * Idm

            en = []
            st = []

            # Computing the levels of momentum k
            if method == :krylov
                en, st, _ = eigsolve(Hk, size(H, 1), nlevels, :SR, ComplexF64; ishermitian = true)
                en = real(en)
            elseif method == :exact
                Hk = Hermitian(Hk)
                en, st = eigen(Hk)
                st = collect(eachcol(st))
            end

            # We initialize the vectors of energies and states at the momentum k/π
            energies[f] = []
            states[f] = []

            for i in 1:nlevels
                if i > length(en) || en[i] > Emax + ΔE/(2 * div_fraction)
                    break
                end 
                push!(energies[f], en[i])
                push!(states[f], st[i])
            end
        end

        # We plot the dispersion relation
        # We extract the couple (momentum, energy) for each state
        k = []
        E = []
        for (kl, en) in energies
            for e in en
                push!(k, kl * π)
                push!(E, e)
            end
        end

        plt = Plots.plot(k, E, seriestype = :scatter, label = "Dispersion relation")
    end

    return args
end

function selectstates(args...)
    # If the :defaultlevel keyword is provided, we select the states at that level
    if haskey(kwargs, :defaultlevel)
        defaultlevel = kwargs[:defaultlevel]
        states = Dict{Rational, Vector{Vector{Complex}}}()
        for (kl, en) in disprel
            states[kl] = [en[defaultlevel]]
        end
        return states
    end
end

function wannier(args...)
    # Def. The (unique) groundstate
    groundstate = states[0 // 1][1]

    # Def. The single-particle dispersion relation, a dictionary momentum => state
    singleparticle = Dict{Rational, Vector{Complex}}()
    
    for (kl, st) in states
        if kl != 0 // 1
            singleparticle[kl] = st[1]
        end
    end
end

# Step 2: l'utente decide di scegliere una specifica particella (una banda), e una volta
# scelta l'algoritmo calcola la Wannier function (sfruttando le simmetrie possibilmente).

# Step 3: l'utente conferma che la Wannier va bene, e da questa calcola l'operatore di
# creazione della particella, lo applica al vuoto DMRG e calcola le funzioni di
# correlazione da cui si possono estrarre le simulazioni di scattering.

let
    # We define the local Hamiltonian
    λ = 0.0

    # We calculate the dispersion relation
    disprel(H0_glueballs(gE = λ, gB = 1 - λ); L0 = 3, Lrange = [9], method = :krylov, nlevels = 20)
end