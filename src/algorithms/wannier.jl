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
export find_wannier