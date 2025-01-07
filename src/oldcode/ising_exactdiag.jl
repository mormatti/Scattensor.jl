# In this snippet, we define the Hamiltonian of the Ising model with transverse and longitudinal fields.

params = Dict(
    :L => 10, # The number of sites of the chain
    :J => 1.0, # The coupling constant of the spins
    :hx => 0.0, # The transverse field
    :hz => 0.0 # The longitudinal field
)

locals = Dict(
    :σx => [0 1; 1 0], # The Pauli matrix σˣ
    :σz => [1 0; 0 -1] # The Pauli matrix σᶻ
)

# We define the Ising Hamiltonian
function hamilt()
    σx, σz = locals[:σx], locals[:σz]
    L, J, hx, hz = params[:L], params[:J], params[:hx], params[:hz]
    𝒪(args...) = matrix_product_local_operators(L, 2, args...)
    return sum(𝒪((J * σz,j),(σz,j+1)) + 𝒪((hx * σx + hz * σz, j)) for j in 1:L)
end


function hamilt(j)
    σx, σz = locals[:σx], locals[:σz]
    L, J, hx, hz = params[:L], params[:J], params[:hx], params[:hz]
    𝒪(args::Pair{Matrix{ComplexF64}, Int}...) = matrix_product_local_operators(L, 2, args...)

    jm = (j-1) ↻ L
    j = j ↻ L
    jp = (j+1) ↻ L

    return 𝒪((J/2 * 𝛔ᶻ,jm),(𝛔ᶻ,j)) + 𝒪((J/2 * 𝛔ᶻ,j),(𝛔ᶻ,jp)) + 𝒪((hˣ * 𝛔ˣ + hᶻ * 𝛔ᶻ,j))
end
