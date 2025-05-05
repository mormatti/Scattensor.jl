# Selecting the groundstate and the first band
bpsi0 = pop_groundstate!(drel)
psi0 = convert_to(MPS, wavefunction(bpsi0), d, Li; cutoff = 1e-12)
E0 = energy(bpsi0)
firstband = get_firstband(drel)

# Computing and Plotting the maximally localized Wannier function
psiw = wannier_symmetric(firstband, H0, L0, Hi, Li, Ti, Ri, d, E0)
veceng = local_exp_value(H0, psiw, L0, d, Li, addconst = -E0/Li)
Plots.plot(log.(abs.(veceng)))