Lbig = 50
sitesL = siteinds(d, Lbig)

# We consider a large system
sitesL = siteinds(d, Lbig)

# The MPO of the local Hamiltonian H0
H0_mpo = mpo_from_matrix(Matrix(H0), d)

# The MPO of the Hamiltonian of the large system
H_mpo = summation_local(H0_mpo, Lbig)
replace_siteinds!(H_mpo, sitesL)

# An initial random MPS to start DMRG
initial_random_mps = random_mps(sitesL)
replace_siteinds!(initial_random_mps, sitesL)

# We compute the groundstate with DMRG
E0_big, psi0_mps_big = dmrg(H_mpo, initial_random_mps, nsweeps = 40)

println("Computing the local basis...")
dim = Lbig - length(creator)
basis = []
for j in 1:dim
    print(" ", j)
    localcreator = insert_local(j - 1, creator, dim - j + 1)
    replace_siteinds!(localcreator, sitesL)
    basisel = apply(localcreator, psi0_mps_big)
    truncate!(basisel, cutoff = 10e-13)
    normalize!(basisel)
    push!(basis, basisel)
end
print(".")

# We compute the correlation matrix and the 
Jstart = (dim + 1) ÷ 2
range = (dim + 1) ÷ 4
CorrMat = []
HMat = []
println("Computing the matrix elements...")
for j in Jstart:(Jstart + range)
    print(" ", j)
    push!(CorrMat, inner(basis[Jstart], basis[j]))
    push!(HMat, inner(basis[Jstart]', H_mpo, basis[j]))
end
print(".")

println("Corr Mat:")
println(CorrMat)
println("HMat")
println(HMat)

Ene(k) = 2 * sum([real(HMat[j+1]) * cos(j * k) for j in 1:6]) + E0 + info["centroid"]
Vel(k) = 2 * sum([real(HMat[j+1]) * (-j) * sin(j * k) for j in 1:6])
Dis(k) = 2 * sum([real(HMat[j+1]) * (-j^2) * cos(j * k) for j in 1:6])

Plots.plot(band)
Plots.plot!(Ene, 0:0.01:π)
# Minimum k dispersion: k \sim 1.712, with speed 0.2333