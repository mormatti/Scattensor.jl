# In this snippet, we define the Hamiltonian of the Ising model with transverse and longitudinal fields.

"""
Generate the Hamiltonian of the Ising model with transverse and longitudinal fields.

Inputs:
- `L` is the number of sites of the chain;
- `J` is the coupling constant of the spins;
- `hˣ` is the transverse field;
- `hᶻ` is the longitudinal field;
- `periodic` is a boolean that indicates if the chain is periodic.
"""
function opsum_ising_hamiltonian(
L::Int64;
J::Float64 = 1.0,
hˣ::Float64 = 0.0,
hᶻ::Float64 = 0.0,
periodic::Bool = false
)

os = OpSum()

# We add the interaction term
if J != 0.0
    for j ∈ 1:L-1
        os += 4 * J, "Sz", j, "Sz", j+1
    end
end

# We add the transverse field term
if hˣ != 0.0
    for j ∈ 1:L
        os += 2 * hˣ, "Sx", j
    end
end

# We add the longitudinal field term
if hᶻ != 0.0
    for j ∈ 1:L
        os += 2 * hᶻ, "Sz", j
    end
end

# We add the periodic boundary conditions
if periodic
    os += 4 * J, "Sz", L, "Sz", 1
end

return os
end

"""
Generate the local Hamiltonian of the Ising model with transverse and longitudinal fields.

Inputs:
- `L` is the number of sites of the chain;
- `j` is the position of the local Hamiltonian;
- `J` is the coupling constant of the spins;
- `hˣ` is the transverse field;
- `hᶻ` is the longitudinal field;
- `periodic` is a boolean that indicates if the chain is periodic.
"""
function opsum_ising_local_hamiltonian(
L::Int64,
j::Int64;
J::Float64 = 1.0,
hˣ::Float64 = 0.0,
hᶻ::Float64 = 0.0,
periodic::Bool = false
)

# We assert that the site is within the chain
@assert 1 ≤ j ≤ L

# We construct the operator sum
os = OpSum()
if j == 1
    if periodic
        os += 0.5 * 4 * J, "Sz", L, "Sz", 1
    end
    os += 0.5 * 4 * J, "Sz", 1, "Sz", 2
    os += 2 * hˣ, "Sx", 1
    os += 2 * hᶻ, "Sz", 1
elseif j == L
    os += 0.5 * 4 * J, "Sz", L-1, "Sz", L
    if periodic
        os += 0.5 * 4 * J, "Sz", L, "Sz", 1
    end
    os += 2 * hˣ, "Sx", L
    os += 2 * hᶻ, "Sz", L
else
    os += 0.5 * 4 * J, "Sz", j-1, "Sz", j
    os += 0.5 * 4 * J, "Sz", j, "Sz", j+1
    os += 2 * hˣ, "Sx", j
    os += 2 * hᶻ, "Sz", j
end
end

λ = 0.2
L = 11
nₑ = 12

# Generation of the Hamiltonian MPO
sites = siteinds("S=1/2", L)
𝐇 = MPO(opsum_ising_hamiltonian(L; J = -λ, hˣ = -(1-λ), periodic = true), sites)

# DMRG parameters
χ₀ = 200 # The initial bond dimension
χₘ = 200 # The maximum bond dimension
sₘ = 10000 # The maximum number of sweeps
ϵᴰ = 10^(-15) # The SVD cutoff
observer = DMRGObserver(energy_tol = 10^(-10)) # The energy tolerance

# We run the DMRG algorithm in order to find the ground state and
# the first nₑ excited states
E₀, ψ₀ = dmrg(
    𝐇, 
    random_mps(sites; linkdims = χ₀); 
    nsweeps = sₘ, 
    maxdim = χₘ, 
    cutoff = ϵᴰ, 
    observer = observer
    )

ψ::Vector{MPS} = []
E = []
for α ∈ 1:nₑ
    Eₙ, ψₙ = dmrg(
        𝐇, 
        [ψ₀; ψ], 
        randomMPS(sites; linkdims = χ₀); 
        nsweeps = sₘ, 
        maxdim = χₘ, 
        cutoff = ϵᴰ, 
        observer = observer)
    push!(E, Eₙ)
    push!(ψ, ψₙ)
end

# We compute the operators in the subspace of found states
𝐇ᵣ = [[inner(ψ1', 𝐇, ψ2) for ψ1 in ψ] for ψ2 in ψ]
𝐓ᵣ = [[inner(ψ1', translate(ψ2)) for ψ1 in ψ] for ψ2 in ψ]
# 𝐑ᵣ = [[inner(ψ1', reflect(ψ2)) for ψ1 in ψ] for ψ2 in ψ]
# We convert the vector of vectors into a matrix of ComplexF64
𝐇ᵣ = Matrix{ComplexF64}(reduce(hcat, 𝐇ᵣ))
𝐓ᵣ = Matrix{ComplexF64}(reduce(hcat, 𝐓ᵣ))
# 𝐑ᵣ = Matrix{ComplexF64}(reduce(hcat, 𝐑ᵣ))
# We compute the simultaneous diagonalization and we plot the dispersion relation
k, ℰ, v = simultaneous_diagonalization_HU(𝐇ᵣ, 𝐓ᵣ)
print("Energies from TN: ", ℰ)
plot_dispersion_relation(k, ℰ, file_name = "disp_TN.png", aspect_ratio = 1.0, points_color = RGBA{Float64}(0.1, 0.1, 1, 1), background_color = RGBA{Float64}(1, 1, 1, 0), marker_size = 1, dpi = 500)

# We compute the exact matrices
𝐇ₑ = matrix_ising_hamiltonian(L, -λ, -(1-λ), 0.0)

𝐓ₑ = translation_operator(2, L)
# We diagonalize 𝐇ₑ and we select the first nₑ states
(Eₑ, 𝛙ₑ) = eigen(𝐇ₑ)
# We write the first energy
𝛙ₑ = [𝛙ₑ[:,i] for i in 1:nₑ]
# We construct the restricted matrices
𝐇ₑᵣ = [[𝛙1' * (𝐇ₑ * 𝛙2) for 𝛙1 in 𝛙ₑ] for 𝛙2 in 𝛙ₑ]
𝐓ₑᵣ = [[𝛙1' * (𝐓ₑ * 𝛙2) for 𝛙1 in 𝛙ₑ] for 𝛙2 in 𝛙ₑ]
# We convert the vector of vectors into a matrix of ComplexF64
𝐇ₑᵣ = Matrix{ComplexF64}(reduce(hcat, 𝐇ₑᵣ))
𝐓ₑᵣ = Matrix{ComplexF64}(reduce(hcat, 𝐓ₑᵣ))
# We compute the simultaneous diagonalization and we plot the dispersion relation
ke, Ee, ψe = simultaneous_diagonalization_HU(𝐇ₑᵣ, 𝐓ₑᵣ)
plot_dispersion_relation(ke, Ee; file_name = "disp_exact.png", aspect_ratio = 1.0, points_color = RGBA{Float64}(1, 0.1, 0.1, 1), background_color = RGBA{Float64}(1, 1, 1, 0), marker_size = 2, dpi = 500)

println("Groundstate energy TN: ", E₀) 
println("Groundstate energy ED: ", Eₑ[1])
println("Error gs energies: ", abs(E₀ - Eₑ[1]) / E₀)