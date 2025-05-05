ϵ = 1e-12

# Creating the MPO wannier operator
psiw = convert_to(MPS, psiw, d, L; cutoff = ϵ)
W = product_outer(psiw, psi0)
truncate!(W, cutoff = ϵ)
Lw = 7
W = partial_trace(W, 2, 8)
truncate!(W, cutoff = ϵ)

L = 100 #100
k0 = π/2
j0 = L/4
sigma = L/100 * 3
χmax = 100 #100
dt = 1
Dt = 240 #400

# We consider a large system
H0 = convert_to(MPO, Matrix(H0), d, 3; cutoff = ϵ)
H = summation_local(H0, d, L; cutoff = ϵ)
psir = random_mps(siteinds(d, L))
substitute_siteinds!(H, psir)
E0, psi0 = dmrg(H, psir, nsweeps = 20)

# We create the wavepacket creator for this large system
f(j) = exp(im * k0 * j) * exp(-(j - (j0 - Lw/2))^2 / (2 * sigma^2))
psi1 = summation_local(W, d, L; func = f, cutoff = ϵ)
truncate!(psi1, cutoff = ϵ, maxdim = χmax)

f(j) = exp(im * -k0 * j) * exp(-(j - (j0 + L/2 - Lw/2))^2 / (2 * sigma^2))
psi2 = summation_local(W, d, L; func = f, cutoff = ϵ)
truncate!(psi2, cutoff = ϵ, maxdim = χmax)

# # We apply the creator
psi = product_matricial(psi1, psi0)
normalize!(psi)
psi = product_matricial(psi2, psi)
normalize!(psi)

println("Done!")

data = Dict()

tdvp_time_evolution!(data, H, psi, dt, Dt, H0, maxdim = χmax, cutoff = ϵ)

println("Computing local expectation values...")
energymatrix = local_expvals(states, H0, d)
Plots.heatmap(energymatrix)