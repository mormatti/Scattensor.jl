using ITensors, ITensorTDVP, ITensorGLMakie
using LinearAlgebra, LinearSolve
using SparseArrays
using Optim
using Plots

function WannierTest()
    w = 4
    λ = 0.8
    L = 4 * w + 1 # The number of sites of the intermediate size chain
    nₑ = L # Number of excited levels to be computed, L if all

    # We create the sites, the Hamiltonian and the list of local Hamiltonians
    sites = siteinds("S=1/2", L)
    𝐇 = MPO(Opsums.ising(L; J = 1-λ, hˣ = λ, periodic = true), sites)
    𝐡 = [MPO(Opsums.ising_local(L, j; J = 1-λ, hˣ = λ, periodic = true), sites) for j ∈ 1:L]

    # # We run the DMRG algorithm in order to find the ground state and the first excited states
    χ₀ = 200 # The initial bond dimension of DMRG
    χₘ = 200 # The maximum bond dimension of DMRG
    sₘ = 10000 # The maximum number of sweeps of DMRG
    ϵᴰ = 10^(-15) # The cutoff of DMRG
    observer = DMRGObserver(energy_tol = 10^(-10))
    E₀, 𝛙₀ = dmrg(𝐇, randomMPS(sites; linkdims = χ₀); nsweeps = sₘ, maxdim = χₘ, cutoff = ϵᴰ, observer = observer)

    # We define important functions
    # energy_density(𝛟,j) = real(inner(𝛟', 𝐡[j], 𝛟) - inner(𝛙₀', 𝐡[j], 𝛙₀))
    parity(𝛟) = inner(𝛟, reflect(𝛟))
    cosk(𝛟) = inner(𝛟, translate(𝛟))

    # We run the DMRG algorithm in order to find the first nₑ excited states
    𝛙::Vector{MPS} = []
    E = []
    for α ∈ 1:nₑ
        Eₙ, 𝛙ₙ = dmrg(𝐇, [𝛙₀; 𝛙], randomMPS(sites; linkdims = χ₀); nsweeps = sₘ, maxdim = χₘ, cutoff = ϵᴰ, observer = observer)
        push!(E, Eₙ)
        push!(𝛙, 𝛙ₙ)
        # plot!([energy_density(𝛙ₙ,j) for j ∈ 1:L], label = "|ψ$α⟩")
    end

    # We compute the symmetrized states and we write the parity
    r = [parity(𝛙[α]) for α ∈ eachindex(𝛙)]
    cosk₀ = cosk(𝛙₀)
    𝛙ₛ = [(𝛙[α] + reflect(𝛙[α]))/(norm(𝛙[α] + reflect(𝛙[α]))) for α ∈ eachindex(𝛙)]
    rₛ = [parity(𝛙ₛ[α]) for α ∈ eachindex(𝛙)]
    println("Parity of the symmetrized states: ")
    println(rₛ)

    # We compute the momenta and we write it
    k₀ = acos(clamp(cosk₀, -1, 1))
    k = [acos(clamp(cosk(𝛙ₛ[α]), -1, 1)) for α ∈ eachindex(𝛙)]
    println("Momenta of the symmetrized states: ", k)

    # We create a scatterplot of the energies vs momenta
    scatter([k₀], [E₀], label = "E0")
    scatter!(k, E, label = "E(k)")
    savefig("dispersion.png")


    # Constructing the matrices
#=  H0 = MPO(opSumIsing(L, λ; n = 0, prefactor = 1/L), sites)
    H1 = MPO(opSumIsing(L, λ; n = 1, prefactor = 1/L), sites)
    H2 = MPO(opSumIsing(L, λ; n = 2, prefactor = 1/L), sites)
    A0 = inner(𝛙₀', H0, 𝛙₀)
    A1 = inner(𝛙₀', H1, 𝛙₀)
    A2 = inner(𝛙₀', H2, 𝛙₀)
    B0 = reduce(hcat, [[inner(𝛙[β]', H0, 𝛙[α]) for α ∈ eachindex(𝛙)] for β ∈ eachindex(𝛙)])
    B1 = reduce(hcat, [[inner(𝛙[β]', H1, 𝛙[α]) for α ∈ eachindex(𝛙)] for β ∈ eachindex(𝛙)])
    B2 = reduce(hcat, [[inner(𝛙[β]', H2, 𝛙[α]) for α ∈ eachindex(𝛙)] for β ∈ eachindex(𝛙)])
 =#
    # Assuming the matrix real, we find the eigenvalues and eigenvectors
    # eig = eigen(B2)
    # vals = eig.values
    # vecs = eig.vectors
    # zꜜ = vecs[:,1]
    # println(vals[1])
    # println("Coefficients: ", zꜜ)

    # Following steps if we cannot assume the matrix real
    #= function ℰ(n,θ)
        A = n==0 ? A0 : n==1 ? A1 : A2
        B = n==0 ? B0 : n==1 ? B1 : B2
        ℑ = eachindex(θ)
        return real(-L * A + sum(exp(im * (θ[α] - θ[β])) * B[α,β] for α ∈ ℑ, β ∈ ℑ))
    end
    σ²(θ) = ℰ(2,θ) / ℰ(0,θ) # - (ℰ(1,θ) / ℰ(0,θ))^2
    θ₀ = zeros(Float64, nₑ)
    result = optimize(σ², θ₀)
    θꜜ = result.minimizer
    println(result) =#

    # We sum all the 𝛙ₙ modulated by the z
    # 𝓌 = sum([zꜜ[α] * 𝛙[α] for α ∈ eachindex(𝛙)], cutoff = 10^(-15))

    # plot!([energy_density(𝓌,j) for j ∈ 1:L], label = "|W⟩")
    
    # plot!([log(abs(energy_density(𝓌,j))) for j ∈ 1:L], label = "|W⟩")

    # savefig("WannierLog.png")

    # We compute the outer product between the Wannier and the groundstate
    # A = outer(𝓌', 𝛙₀)
end