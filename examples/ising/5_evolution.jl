# Select the k0 which guarantees v > 0 and minimize the dispersion
k0 = 1.712

# Parameters obtained from k0
lw = lws[1]
v0 = Vel(k0)
d0 = Dis(k0)
σ0 = Lbig / 36
Dt = Lbig / v0
j0 = 2 * Lbig / 6
χmax = 100 #100
dt = Dt / 200

println("Generating the left- creator...")
# We create the wavepacket creator for this large system
wavef(j,jbar,kbar,σbar) = exp(im * kbar * j) * exp(-(j - (jbar - lw))^2 / (2 * σbar^2))
creator_left = summation_local(creator, Lbig; convolution = (j -> wavef(j, j0 - Lbig/3, k0, σ0)))
creator_right = summation_local(creator, Lbig; convolution = (j -> wavef(j, j0 + Lbig/3, -k0, σ0)))

println("Initializing the state...")
substitute_siteinds!(psi0_mps_big, creator_left)
psi = normalize(apply(creator_left, psi0_mps_big))
substitute_siteinds!(psi0_mps_big, creator_right)
psi = normalize(apply(creator_right, psi0_mps_big))

println("Computing time evolution...")
data = Dict()
tdvp_time_evolution!(data, H_mpo, psi, dt, Dt, H0_mpo, maxdim = χmax, cutoff = ϵcutoff)

println("Computing local expectation values...")
energymatrix = local_expvals(states, H0, d)
Plots.heatmap(energymatrix)