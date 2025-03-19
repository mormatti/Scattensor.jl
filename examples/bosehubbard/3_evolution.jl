ϵ = 1e-12

# Creating the MPO wannier operator
psiw = convert_to(MPS, psiw, d, L; cutoff = ϵ)
W = product_outer(psiw, psi0)
truncate!(W, cutoff = ϵ)
W = partial_trace(W, 2, 12)
Lw = 10
truncate!(W, cutoff = ϵ)

L = 100 #100
k0 = π/2
j0 = L/4
sigma = L/100 * 2
χmax = 100 #100
dt = 0.5
Dt = 400 #400

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

# We apply the creator
psi = product_matricial(psi1, psi0)
normalize!(psi)
psi = product_matricial(psi2, psi)
normalize!(psi)

# We compute the time evolution and we observe local operators
states, info = tdvp_time_evolution(H, psi, dt, Dt, maxdim = χmax, cutoff = ϵ)
println("Simulation ended. Total size of states: $(Base.summarysize(states) / 1000000) MB ")

println("Computing local expectation values...")
energymatrix = local_expvals(states, H0, d)
Plots.heatmap(energymatrix)