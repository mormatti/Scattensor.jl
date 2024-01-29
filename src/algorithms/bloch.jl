"""
Symultaneously diagonalize an Hermitian matrix 𝐇 and a unitary matrix 𝐔 such that
they commute, i.e. [𝐇,𝐔] = 0.

Assumptions:
- `𝐇` must be translationally invariant, i.e. [𝐇,𝐔] = 0.

Inputs:
- `𝐇` is an Hermitian matrix, i.e. 𝐇 = 𝐇†;
- `𝐔` is a unitary matrix, i.e. 𝐔𝐔† = 𝐔†𝐔 = 𝟙.

Outputs:
- `u` is the `Vector` of the phases (angles) of the eigenvalues of 𝐔;
- `h` is the `Vector` of real eigenstates associated to 𝐇;
- `𝛙` is the `Vector` of (common) eigenvectors associated to 𝐇 and 𝐔.
"""
function bloch_states(
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
export bloch_states


"""
Plots the dispersion relation of the system.

Inputs:
- `k` is the `Vector` of momenta of the eigenstates;
- `ℰ` is the list of energies of the eigenstates;
- `filename` is the name of the file (with extension) where the plot will be saved;
- `aspectRatio` is the ratio between the width and the height of the plot (default 1:1).
"""
function plot_dispersion_relation(
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
export plot_dispersion_relation

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
export firstBandStates

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
export energy_density_array