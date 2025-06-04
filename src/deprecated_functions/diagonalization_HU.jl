# Useful but old function...
""" Simultaneously diagonalize an Hermitian matrix 𝐇 and a unitary matrix 𝐔 such that they commute, i.e. [𝐇,𝐔] = 0.

    ## Assumptions
    - `𝐇` must be translationally invariant, i.e. [𝐇,𝐔] = 0.

    ## Inputs
    - `𝐇` is an Hermitian matrix, i.e. 𝐇 = 𝐇†;
    - `𝐔` is a unitary matrix, i.e. 𝐔𝐔† = 𝐔†𝐔 = 𝟙.

    ## Outputs
    - `u` is the `Vector` of the phases (angles) of the eigenvalues of 𝐔;
    - `h` is the `Vector` of real eigenstates associated to 𝐇;
    - `𝛙` is the `Vector` of (common) eigenvectors associated to 𝐇 and 𝐔.
    """
function simultaneous_diagonalization_HU(
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
        @assert ishermitian(𝐇) "𝐇 is not hermitian. |𝐇† - 𝐇|/|𝐇| = $(norm(𝐇' - 𝐇)/norm(𝐇))."
    end 

    if check_unitarity
        @assert 𝐔 * 𝐔' ≈ 𝐔' * 𝐔 ≈ I "𝐔 is not unitary. |𝐔†𝐔 - I| = $(norm(𝐔' * 𝐔 - I))."
    end

    if check_translational_invariance
        @assert 𝐇 * 𝐔 ≈ 𝐔 * 𝐇 "𝐇 is not translational invariant. |𝐇𝐔 - 𝐔𝐇|/|𝐇𝐔| = $(norm(𝐇𝐔 - 𝐔𝐇)/norm(𝐇𝐔))."
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