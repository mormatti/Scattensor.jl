Lbig = 20
χmax = 50

# We consider a large system
H0mpo = convert_to(MPO, Matrix(H0), d, 3; cutoff = ϵcutoff)
Hmpo = summation_local(H0mpo, d, Lbig; cutoff = ϵcutoff)
psir = random_mps(siteinds(d, Lbig))
substitute_siteinds!(Hmpo, psir)
E0big, psi0big = dmrg(Hmpo, psir, nsweeps = 20)