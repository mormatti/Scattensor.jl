using Scattensor
using Revise

# We select the groundstate and the first band from the dispersion relation
band = get_firstband(disprel)
# band2 = pop_firstband!(disprel)
# band3 = pop_firstband!(disprel) # Just to check that we have more than 2 bands
# bands = [band1] # Just to check that we have more than 2 bands

# We define the local hamiltonian in the central site
Llat = Int((Li - L0)/2) # The number of sites to add on each side of the local hamiltonian
# Hc = kron_power(id, Llat) ⊗ H0 ⊗ kron_power(id, Llat)
Hc_mpo = insert_local(Llat, H0_mpo, Llat)
replace_siteinds!(Hc_mpo, sitesl)

# wvec, info = wannier_symmetric(band, E0, Hc, Pi)
wmps, info = wannier_symmetric(band, E0, x -> apply(Hc_mpo, x), x -> apply_reflection(x), (x,y) -> inner(x,y))
truncate!(wmps, cutoff = 1e-13)

# wmps = mps_from_vector(wvec, d) # The wannier in MPS formi
# replace_siteinds!(wmps, sitesl)

Plots.plot(Dict(i => log10(abs(info["density"][i])) for i in keys(info["density"])))