ϵcutoff = 1e-10

# We create the wavepacket creator for this large system
psiw = convert_to(MPS, psiw, d, Li; cutoff = ϵcutoff)
W = product_outer(psiw, psi0)
truncate!(W, cutoff = ϵcutoff)

W = partial_trace(W, 5, 7)
truncate!(W, cutoff = ϵcutoff)