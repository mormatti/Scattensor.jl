using LinearAlgebra
using ITensors

lw = 5 # Set this
jc = Int((Li - 1)/2 + 1) # The central site

# This works!

# A = partial_trace(outer(ωmps', wmps), jc - lw, jc + lw)
# Svd = svd(matrix_from_mpo(A))
# println("Nuclear norm = ", sum(Svd.S))
# println("Singular values: ", Svd.S[1:3])
# creator = mpo_from_matrix(Svd.V * Svd.U', d)
# creator = mpo_from_matrix((Svd.V[:,1]) * (Svd.U')[:,1]', d)

# But magically also this! (When the wannier is nice)
# creator = partial_trace(outer(wmps', ωmps), jc - lw, jc + lw)
# truncate!(creator, cutoff = 1e-10)

Trwω = partial_trace(truncate(outer(wmps', ωmps), cutoff = 1e-10), jc - lw, jc + lw)
sitesmall = siteinds_main(Trwω)
Trωw = partial_trace(truncate(outer(ωmps', wmps), cutoff = 1e-10), jc - lw, jc + lw)
replace_siteinds!(Trωw, sitesmall)
Heff1 = truncate(-product(Trwω, Trωw), cutoff = 1e-10)
Heff2 = truncate(-product(Trωw, Trwω), cutoff = 1e-10)
μ1, v1 = dmrg(Heff1, random_mps(sitesmall, linkdims = 10), nsweeps = 20, maxdim = 500, cutoff = 1e-10)
μ2, v2 = dmrg(Heff2, random_mps(sitesmall, linkdims = 10), nsweeps = 20, maxdim = 500, cutoff = 1e-10)
print("Sigma max = $(√(-μ0))")
creator = truncate(outer(v1', v2))

wlat = Int((Li - 2 * lw - 1)/2)
creator_on_l = insert_local(wlat, creator, wlat)
replace_siteinds!(creator_on_l, sitesmall)
inner(wmps, apply(creator_on_l, ωmps))