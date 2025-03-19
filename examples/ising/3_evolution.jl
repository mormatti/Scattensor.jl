# Creating the MPO wannier operator
psiw = convert_to(MPS, psiw, d, L; cutoff = 1e-12)
W = product_outer(psiw, psi0)
truncate!(W, cutoff = 1e-12)
W = partial_trace(W, 3, 9)
truncate!(W, cutoff = 1e-12)

# We consider a large system
L = 30
H0 = convert_to(MPO, Matrix(H0), d, 3; cutoff = 1e-12)
H = summation_local(H0, d, L; cutoff = 1e-12)
psir = random_mps(siteinds(d, L))
substitute_siteinds!(H, psir)
E0, psi0 = dmrg(H, psir, nsweeps = 10)

# We create the wavepacket creator for this large system
kbar = π/2
jbar = 10
sigma = 20
f(j) = exp(im * kbar * j) * exp(-(j - jbar)^2 / (2 * sigma))
psi = summation_local(W, d, L; func = f, cutoff = 1e-12)
truncate!(psi, cutoff = 1e-12, maxdim = 30)

# We apply the creator
psi = product_matricial(psi, psi0)
normalize!(psi)

# We compute the time evolution and we observe local operators
dt = 0.1
Dt = 3
states, info = tdvp_time_evolution(H, psi, dt, Dt, maxdim = 30, cutoff = 1e-12)

energymatrix = local_expvals(states, H0, d)
Plots.heatmap(energymatrix)