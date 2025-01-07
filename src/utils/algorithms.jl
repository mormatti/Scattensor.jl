"""
    Finds the bloch states of the first band.

    Inputs:
    - `k` is the `Vector` of momenta of the eigenstates;
    - `ℰ` is the list of energies of the eigenstates;
    - `𝛙` is the `Matrix` of eigenvectors of the eigenstates.

    Assumptions:
    - The groundstate must be non-degenerate.
    """
function first_band_states(
    k::Vector{Float64},
    ℰ::Vector{Float64},
    𝛙::Vector{Vector{ComplexF64}},
    L::Integer
    )::Tuple{Vector{Float64}, Vector{Float64}, Vector{Vector{ComplexF64}}}

    @debug "Computing the first band states..."

    # Constructing the vector of tuples (k, ℰ, 𝛙)
    𝓣 = [(k[i], ℰ[i], 𝛙[i]) for i in eachindex(k)]

    # Selecting the groundstate and deleting it
    j = argmin([t[2] for t in 𝓣]) # the groundstate index for 𝓣
    deleteat!(𝓣, j) # Deleting the element of states at that index

    # Selecting the first band states
    # The list of the tuple states of the first band
    𝓣′ = []
    while length(𝓣) > 0
        # The index and momentum of the state with the minimum energy
        j′ = argmin([t[2] for t in 𝓣])
        k′ = 𝓣[j′][1]
        # Filling 𝓣′ with the state with the minimum energy
        push!(𝓣′, 𝓣[j′])
        # deleting all the elements with same (similar) momentum
        𝓣 = [t for t in 𝓣 if abs(k′ - t[1]) > π/L]
        # we repeat the process until 𝓣 is empty
    end

    # Transforming the a vector of tuples 𝓣′ into a tuple of vectors 𝓑′
    # the vector of momenta of the first band
    k = [t[1] for t in 𝓣′]
    # the vector of energies of the first band
    ℰ = [t[2] for t in 𝓣′]
    # the matrix of eigenvectors of the first band
    𝛙 = [t[3] for t in 𝓣′]

    # Returning the tuple of bloch states of the first band
    return k, ℰ, 𝛙
end

""" 
    Plots the energy density profile of the state `𝛙`.

    Inputs:
    - `𝛙` is a generic qualtum state;
    - `𝐡₁` is the local Hamiltonian of the first site;
    - `𝐓` is the translation operator of the system;
    - `L` is the number of sites of the chain.

    Optional inputs:
    - `ℰ₀` is the groundstate energy of the system (default 0.0).

    Outputs:
    - The `Vector{Real}` of the energy density profile of the state `𝛙`.
    """
function energy_density_array(
    𝛙::Vector{ComplexF64}, 
    𝐡₁::Matrix{ComplexF64}, 
    𝐓::Matrix{ComplexF64}, 
    L::Integer;
    ℰ₀::Float64 = 0.0,
    )::Vector{Float64}

    @debug "Computing the energy density profile..."

    # We initialize the array
    arr::Vector{Float64} = [real(𝛙' * 𝐡₁ * 𝛙 - ℰ₀/L)]

    # The loop to compute all the others energies
    for _ in 2:L
        𝛙 = 𝐓' * 𝛙
        push!(arr, real(𝛙' * 𝐡₁ * 𝛙 - ℰ₀/L))
    end

    return arr
end

"""
    Interpolates a state 𝛙₁ from a state 𝛙₀, using a linear combinations of compositions
    of local operators 𝐋ⱼ.
    """
function interpolate(
    𝛙₁::Vector{ComplexF64}, # The vector to interpolate
    𝛙₀::Vector{ComplexF64}, # The vector from which 𝛙₁ is interpolated
    𝐓::Matrix{ComplexF64}, # The translation operator of the system
    j₁::Int64, # The position of the left border of the interpolation range
    j₂::Int64, # The position of the right border of the interpolation range
    𝐋::Vector{Matrix{ComplexF64}}, # A list of d local operators
    d::Int64, # The dimension of the local space
    L::Int64, # The number of sites of the chain
    )
    # We define the identity
    Ide = Matrix{ComplexF64}(I, d, d)

    function generate_𝐋₁(n)
        𝐀 = 𝐋[n]
        for _ in 2:L
            𝐀 = 𝐀 ⊗ Ide
        end
        return 𝐀
    end

    𝐋₁ = [generate_𝐋₁(n) for n in 1:d]

    # Apply 𝐋ⱼ to the state 𝛙 in the most efficient way (translating the state instead 
    # of the operator)
    function apply(n, j, 𝛙)
        𝛟 = deepcopy(𝛙)
        for _ in 1:j
            𝛟 = 𝐓' * 𝛟
        end
        𝛟 = 𝐋₁[n+1] * 𝛟
        for _ in 1:j
            𝛟 = 𝐓 * 𝛟
        end
        return 𝛟
    end

    function apply_to_array(n, j, arr)
        arr′ = deepcopy(arr)
        arr′[j] = n
        return arr′
    end

    𝚿 = [𝛙₀]
    A = [zeros(L)]

    # Find all the possible combinations of the operators starting from 𝛙₀
    # We also find an array representing these applications
    for j in j₁:j₂
        𝚿 = [apply(n,j,𝛙) for n in 0:(d-1) for 𝛙 in 𝚿]
        A = [apply_to_array(n,j,a) for n in 0:(d-1) for a in A]
    end

    # We compute the Hessian matrix
    Hs = [𝚿[α]' * 𝚿[β] for α in eachindex(𝚿), β in eachindex(𝚿)]
    b = [𝚿[α]' * 𝛙₁ for α in eachindex(𝚿)]
    println("Dimension of the Hessian matrix: ", size(Hs))
    println("Rank of the Hessian matrix: ", rank(Hs))

    # We compute the coefficients, calculating the pseudo-inverse of 𝚿
    println("Computing the inverse of the Hessian matrix...")
    c = pinv(Hs) * b

    # We compute the infidelity
    iFi = 1 + sum(c[α]' * c[β] * Hs[α,β] for α in eachindex(c), β in eachindex(c))
                 - sum((c[α]' * b[α] + c[α] * b[α]') for α in eachindex(c))
    println("Infidelity: ", iFi)

    # We create a list of the tuples (coefficients, A)
    listt = [(c[α], A[α]) for α in eachindex(c)]

    # We sort the list by the coefficients
    sort!(listt, by = x->abs(x[1]), rev = true)

    # We print the first 20 elements of listt in a nice way
    # println()
    # println("The first 20 elements are:")
    # println()
    # for α in 1:20
        # println(listt[α][2], " ", abs(listt[α][1]))
    # end

    # we create the array of the abs of the coefficients sorted
    xax = [i for i in 1:length(listt)]
    cp = [abs(x[1]) for x in listt]

    plot(xax, cp, xaxis=:log, yaxis=:log)
    savefig("data/coefficients.png")
end

"""
Assumptions: PBC with translation operator like T^L = I.
Here the ordinary bloch theorem is valid.
"""
function wavenumber(z, L)::Int
    return Int(round(L/2π * angle(z)))
end

# TODO: use Diagonal instead of dictionary
""" Compute the maximally localized Wannier function.

    Let us consider a 1D lattice system with L sites and unit spacing a.
    The Hamiltonian H is given by ∑ⱼhⱼ, where hⱼ is the local Hamiltonian.
    Assuming translational invariance, each hⱼ is the same.

    # Assumptions (temporary)
    - H = ∑ⱼhⱼ (interaction locality)
    - [T,H] = 0 (translational invariance, should be always assumed)
    - [R,H] = 0 (reflection invariance, can be generalized)
    - T^L = I (ordinary boch theorem, no SSB)

    # Arguments
    - `ϵ::Dict` is a sorted dictionary of the form wavenumber => energy. The wavenumber n
    is computed as aLk/2π = φ/a, where φ is the momentum phase of each eigenvector.
    - `h::Matrix` is the local Hamiltonian in the basis of the L bloch states. It is
    a L×L matrix.

    # Definitions
    - The central site: (L+1)/2 for L odd, L/2 for L even

    # Returns

    """
function compute_localized_wannier_coefficients(
    ϵ::Dict, # The sorted dictionary wavenumber => energy (wavenumber = int for T^n = I)
    h::Matrix; # Local hamiltonian in the central site, in the basis of the L bloch states
    E₀::Float64 = 0.0, # The groundstate energy of the system, default 0.0
    reflectinvar::Bool = false # If the system is reflection invariant, default false
    )::Vector{Complex}

    # Define the parameters and basic functions
    L = length(ϵ) # The number of vectors
    @assert size(h) == (L, L) # The local Hamiltonian must be a square matrix of size L×L
    c = (L % 2 == 0) ? (L/2) : (L+1)/2 # The central site
    
    # We shift the energy levels to have the groundstate energy at 0
    ϵ = Dict(k => v - groundstenergy for (k,v) in ϵ)
    h -= groundstenergy/L * I
    
    ϵmean = sum(values(ϵ)) / N # The average of the energy states
    χ(j) = L/π * sin(π/L * (j - c)) # We define the metric

    # We compute the matrix A₂
    A₂ = [[sum([h[α,β] * χ(j)^2 * t[α]^(j-c) * t[β]^(c-j) for j in 1:L]) for α ∈ 1:N] for β ∈ 1:N]
    z(θ) = [exp(i * θα) for θα ∈ θ] # We define a function that given a vector of reals return the vector of phases
    f(θ) = real(z[θ]' * A₂ * z[θ]) # We define the functional

    # We minimize the functional
    result = optimize(f, ones(N-1), NelderMead())
    θₘ = result.minimizer

    wₘ = z(θₘ) / √L # The maximally localized Wannier function
    σₘ = f(θₘ) / ϵmean # The minimum spread of the Wannier function
    Eⱼ = [sum([sum([w[α]' * h[α,β] * t[α]^(j-c) * t[β]^(c-j) * w[β] for α ∈ 1:L]) for β ∈ 1:L]) for j ∈ 1:L]

    # We compute the maximally localized state
    return wₘ, σₘ, Eⱼ
end


""" Most general possible implementation
    - If L is odd, the central site is j = (L+1)/2
    - If L is even, the central site is chosen to be j = L/2
"""
function find_wannier(
    L::Int, # The number of sites of the chain
    ψ::Matrix, # The set of n states (an D×n matrix) which we want to localize
    T::Matrix, # The D×D translation operator of the system
    hc # The local Hamiltonian (a D×D matrix acting on ψ) of the central site ((L+1)/2 if odd, L/2 if even)
    )::Vector{ComplexF64}
    
    @debug "Computing the localized Wannier state..."
    N = length(𝛙) # We define the number of states

    # For each vector in 𝛙 we compute the translated vector by 1 unit
    Tψ = T * ψ

    # For each state in 𝚿 we compute the multiplied by the local Hamiltonian
    ψh = Tψ

    # Computing the overlap matrix Λ between States and StatesH
    Λ = [𝚿[j][α]' * 𝐡𝚿[j][β] for α in 1:N, β in 1:N, j in 1:L]

    # We define the distance function
    Id(j,j′) = 1
    χ(j,j′) = ((Int64(j - j′ - (L-1)/2) ↻ Int64(L)) - (L-1)/2 - 1)
    χ²(j,j′) = χ(j,j′)^2

    # We define the function A₀, A₁ and A₂
    function A(f::Function, θ::Vector{Float64}, j′::Int64)::ComplexF64
        S1 = sum(exp(im * (θ[β] - θ[α])) * Λ[α,β,j] * f(j,j′) for j in 1:L, α in 1:N, β in 1:N)
        S2 = sum(ℰ₀ * f(j,j′) for j in 1:L)
        return (S1 - S2) / L
    end

    # We define auxiliary functions
    A₀(θ, j′) = A(Id, θ, j′)
    A₁(θ, j′) = A(χ, θ, j′)
    A₂(θ, j′) = A(χ², θ, j′)

    # We define the functional to minimize
    function 𝑓(θ::Vector{Float64}, j′::Int64)::Float64
        return real((A₂(θ, j′) / A₀(θ, j′)) - (A₁(θ, j′) / A₀(θ, j′))^2)
    end

    # We define the restricted functional, defined putting the last θ argument to 0
    function 𝑓ᵣ(θr::Vector{Float64})::Float64
        θ = deepcopy(θr)
        push!(θ, 0.0)
        return 𝑓(θ, jᶜ)
    end

    # We minimize the functional
    result = optimize(𝑓ᵣ, ones(N-1), NelderMead())
    θ = result.minimizer
    push!(θ, 0.0)

    println("The minimum is: ", 𝑓(θ, jᶜ))

    # We compute the maximally localized state
    𝓌 = 1/√L * sum(exp(im * θ[j]) * 𝛙[j] for j in 1:N)

    # We compute the maximally localized state
    return 𝓌
end