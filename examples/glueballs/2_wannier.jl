# Selecting the groundstate and the first band
bpsi0 = pop_groundstate!(drel)
psi0 = convert_to(MPS, wavefunction(bpsi0), d, L; cutoff = 1e-12)
E0 = energy(bpsi0)
firstband = get_firstband(drel)

# Computing and Plotting the maximally localized Wannier function
psiw = wannier_symmetric(firstband, H0, L0, H, L, T, R, d, E0)
Plots.plot(log.(abs.(local_exp_value(H0, psiw, L0, d, L, addconst = -E0/L))))

# To plot not in Log:
# Plots.plot(local_exp_value(H0, psiw, L0, d, L, addconst = -E0/L))