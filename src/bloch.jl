using Scattensor
using LinearAlgebra
using Plots
using Optim

"""
Symultaneously diagonalize an Hermitian matrix 𝐇 and a unitary matrix 𝐔 such that
they commute, i.e. [𝐇,𝐔] = 0.

Inputs:
- `𝐇` is an Hermitian matrix, i.e. 𝐇 = 𝐇†;
- `𝐔` is a unitary matrix, i.e. 𝐔𝐔† = 𝐔†𝐔 = 𝟙.

Outputs:
- `u` is the `Vector` of the phases (angles) of the eigenvalues of 𝐔;
- `h` is the `Vector` of real eigenstates associated to 𝐇;
- `𝛙` is the `Vector` of (common) eigenvectors associated to 𝐇 and 𝐔.

Assumptions:
- `𝐇` must be translationally invariant, i.e. [𝐇,𝐔] = 0.
"""
function blochStates(
    𝐇::Matrix{ComplexF64},
    𝐔::Matrix{ComplexF64};
    check_hermiticity::Bool = false,
    check_unitarity::Bool = false,
    check_translational_invariance::Bool = false
    )::Tuple{Vector{Float64}, Vector{Float64}, Vector{Vector{ComplexF64}}}

    @debug "Computing the bloch states from exact diagonalization..."

    # We assert the two matrices have the same size and are square
    @assert size(𝐇)[1] == size(𝐇)[2] == size(𝐔)[1] == size(𝐔)[2] "𝐇 and 𝐔 must be square and having the same size."

    N = size(𝐇)[1] # The dimension of the Hilbert space

    if check_hermiticity
        @assert ishermitian(𝐇) "𝐇 is not hermitian."
    end 

    if check_unitarity
        @assert 𝐔 * 𝐔' ≈ 𝐔' * 𝐔 ≈ I "𝐔 is not unitary."
    end

    if check_translational_invariance
        @assert 𝐇 * 𝐔 ≈ 𝐔 * 𝐇 "𝐇 is not translational invariant."
    end

    # We compute the groundstate energy E₀ and the matrix product 𝐇𝐓,
    # where 𝐇 is shifted by E0 - 1
    h₀ = eigen(𝐇, permute = true).values[1]
    𝐇′𝐔 = (𝐇 - h₀ * I + I) * 𝐔

    # We compute all the eigenvectors and eigenvalues of HT
    (𝜆, 𝛙) = eigen(𝐇′𝐔)
    Nᵥ = length(𝜆) # The number of eigenvalues

    # Extract the phases and moduli from the eigenvalues 𝜆
    (u, h) = (zeros(Nᵥ), zeros(Nᵥ))
    for i in eachindex(𝜆)
        (u[i], h[i]) = (angle(𝜆[i]), abs(𝜆[i]) + h₀ - 1)
    end

    𝛙 = [𝛙[:,i] for i in 1:N]

    return u, h, 𝛙
end


"""
Plots the dispersion relation of the system.

Inputs:
- `k` is the `Vector` of momenta of the eigenstates;
- `ℰ` is the list of energies of the eigenstates;
- `filename` is the name of the file (with extension) where the plot will be saved;
- `aspectRatio` is the ratio between the width and the height of the plot (default 1:1).
"""
function plotDispersionRelation(
    k::Vector{Float64},
    ℰ::Vector{Float64},
    filename::String = "plot.png", 
    aspectRatio::Float64 = 1.0)

    @debug "Plotting the dispersion relation..."

    # Plotting the dispersion relation
    plot(k, ℰ, seriestype=:scatter, markersize=4, legend=false, xlabel="k", ylabel="ℰ")
    plot!(size=(170,170 * aspectRatio), dpi=300) # Setting the size of the plot
    savefig(filename) # Saving the plot in a file
end


"""
Finds the bloch states of the first band.

Inputs:
- `k` is the `Vector` of momenta of the eigenstates;
- `ℰ` is the list of energies of the eigenstates;
- `𝛙` is the `Matrix` of eigenvectors of the eigenstates.

Assumptions:
- The groundstate must be non-degenerate.
"""
function firstBandStates(
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


""" Plots the energy density profile of the state `𝛙`.

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
function energyDensityArray(
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
Computes the maximally localized from a generic set of states.

Inputs:
- `𝛙` is the set of states (a `Vector` of `Vector`s);
- `𝐓` is the translation operator of the system;
- `𝐡` is the local Hamiltonian of the first site;
- `jᶜ` is the localization position of the maximally localized state;
- `ℰ₀` is the groundstate energy of the system;
- `L` is the number of sites of the chain.

Outputs:
- The maximally localized state.
"""
function find_wannier(
    𝛙::Vector{Vector{ComplexF64}},
    𝐓::Matrix{ComplexF64},
    𝐡₁::Matrix{ComplexF64},
    jᶜ::Int64,
    ℰ₀::Float64,
    L::Int64
    )::Vector{ComplexF64}
    
    @debug "Computing the localized Wannier state..."

    # We define the number of states
    N = length(𝛙)

    # For each vector in 𝛙 we compute the translated vector by 1 unit
    𝚿 = [deepcopy(𝛙)]
    𝛙′ = deepcopy(𝛙)
    for _ in 1:(L-1)
        𝛙′ = [𝐓' * 𝜓 for 𝜓 in 𝛙′]
        push!(𝚿, deepcopy(𝛙′))
    end

    # For each state in 𝚿 we compute the multiplied by the local Hamiltonian
    𝐡𝚿 = [deepcopy([𝐡₁ * 𝜙 for 𝜙 in 𝛟]) for 𝛟 in 𝚿]

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


# function find_localized(
#     𝛙::Vector{Vector{ComplexF64}},
#     𝐓::Matrix{ComplexF64},
#     𝐡₁::Matrix{ComplexF64},
#     jᶜ::Int64,
#     L::Int64
#     )::Vector{ComplexF64}
    
#     @debug "Computing the localized Wannier state..."

#     # We define the number of states
#     N = length(𝛙)

#     # For each vector in 𝛙 we compute the translated vector by 1 unit
#     𝚿 = [deepcopy(𝛙)]
#     𝛙′ = deepcopy(𝛙)
#     for _ in 1:(L-1)
#         𝛙′ = [𝐓 * 𝜓 for 𝜓 in 𝛙′]
#         push!(𝚿, deepcopy(𝛙′))
#     end

#     # For each state in 𝚿 we compute the multiplied by the local Hamiltonian
#     𝐡𝚿 = [deepcopy([𝐡₁ * 𝜙 for 𝜙 in 𝛟]) for 𝛟 in 𝚿]

#     # We define the distance function
#     χ(j) = ((Int64(j - jᶜ - (L-1)/2) ↻ Int64(L)) - (L-1)/2 - 1)
#     χ²(j) = χ(j)^2

#     # Computing the overlap matrix Λ between States and StatesH
#     Λ = [sum(χ(j)^2 * 𝚿[j][α]' * 𝐡𝚿[j][β] for j in 1:L) for α in 1:N, β in 1:N]

#     # We select the first (the smallest) eigenvector z of Λ
#     eig = eigen(Λ)
#     z = eig.vectors[:,1]

#     # We compute the maximally localized state
#     lf = sum(z[α] * 𝛙[α] for α in 1:N)
#     lf = lf / norm(lf)

#     # We compute the maximally localized state
#     return lf
# end