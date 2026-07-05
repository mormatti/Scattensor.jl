# MPS spectroscopy port: matrix ↔ MPS agreement + fast OBC pipeline check.
using SparseArrays, LinearAlgebra, Random
using ITensors, ITensorMPS

# kron-convention converter (matrix side: site 1 slowest, site L fastest).
# ITensor arrays are column-major (first index fastest), so site n maps to
# reshaped dimension L−n+1.
function mps_from_kron_vector(v::Vector, sites)
    L = length(sites); d = dim(sites[1])
    T = ITensor(reshape(Vector{ComplexF64}(v), ntuple(_ -> d, L)),
                reverse(collect(sites))...)
    MPS(T, sites; cutoff = 1e-14, maxdim = 4096)
end

@testset "spectroscopy mps" begin
    Random.seed!(11)
    d, L = 2, 8
    σx = sparse([0.0 1.0; 1.0 0.0]); σz = sparse([1.0 0.0; 0.0 -1.0])
    id2 = sparse(1.0I, 2, 2)
    J, h = 1.0, 4.0
    H0 = -J * kron(σz, σz) - 0.5h * (kron(σx, id2) + kron(id2, σx))

    # ---- 1. transition operator: MPS contraction ≡ matrix reshape/QR ----
    sites = siteinds(d, L)
    v = normalize!(randn(ComplexF64, d^L))
    w = normalize!(randn(ComplexF64, d^L))
    ψm = mps_from_kron_vector(v, sites)
    Ωm = mps_from_kron_vector(w, sites)
    for (site, ℓ) in ((3, 3), (1, 2), (5, 4))
        Amat = transition_operator(v, w, ℓ; site)
        Amps = mps_transition_operator(ψm, Ωm, site:(site + ℓ - 1))
        @test norm(Amat - Amps) < 1e-9
    end

    # ---- 2. window application: (A at site) ∘ MPS ≡ embed_operator ∘ vector
    A = randn(ComplexF64, d^3, d^3)
    site = 4
    uvec = embed_operator(A, d, L, site) * w
    umps = mps_apply_window(A, Ωm, site; cutoff = 1e-13, maxdim = 4096)
    ucheck = mps_from_kron_vector(uvec, sites)
    overlap = inner(ucheck, umps) / (norm(ucheck) * norm(umps))
    @test abs(overlap - 1) < 1e-9
    @test abs(norm(umps) - norm(uvec)) < 1e-9

    # ---- 3. deformed trap MPO ≡ matrix deformed Hamiltonian (OBC) ----
    Hm = deformed_hamiltonian(H0, d, L; weight = parabolic_weight(L), pbc = false)
    Hmpo = deformed_hamiltonian_mpo(H0, d, L; weight = parabolic_weight(L))
    ψt = mps_from_kron_vector(v, siteinds_main(Hmpo))
    ev_mpo = inner(ψt', Hmpo, ψt)
    ev_mat = dot(v, Hm * v)
    @test abs(ev_mpo - ev_mat) < 1e-8 * max(1, abs(ev_mat))

    # ---- 4. fast OBC pipeline at L = 12 (deep paramagnet, 1 creator) ----
    L2 = 12
    Hfull = summation_local(mpo_from_matrix(Matrix(H0), d), L2)
    sites2 = siteinds_main(Hfull)
    E0dm, Ω2 = dmrg(Hfull, random_mps(sites2; linkdims = 8);
                    nsweeps = 10, maxdim = [10, 20, 40, 80], cutoff = 1e-11,
                    outputlevel = 0)
    trap = deformed_hamiltonian_mpo(H0, d, L2; weight = parabolic_weight(L2))
    replace_siteinds!(trap, sites2)
    _, seeds = dmrg_seeds(trap, sites2, 2; nsweeps = 10, outputlevel = 0)
    seed = seeds[2]
    a = inner(Ω2, seed)
    seedP = add(seed, -a * Ω2; cutoff = 1e-12)
    normalize!(seedP)
    ℓφ = 3
    tns = [sum(svdvals(mps_transition_operator(seedP, Ω2, s:(s + ℓφ - 1))))
           for s in 2:(L2 - ℓφ)]
    sbest = (2:(L2 - ℓφ))[argmax(tns)]
    A2 = mps_transition_operator(seedP, Ω2, sbest:(sbest + ℓφ - 1))
    A2 ./= norm(A2)
    S, Hc = mps_bulk_couplings(Hfull, Ω2, [A2], L2 ÷ 2, 3)
    ks, bands = wannier_bands(S, Hc; nk = 60, gram_tol = 1e-3)
    eps_exact(k) = 2 * sqrt(J^2 + h^2 - 2J * h * cos(k))
    rel = Float64[]
    for (i, k) in enumerate(ks)
        isnan(bands[i, 1]) && continue
        push!(rel, abs(bands[i, 1] - eps_exact(k)) / eps_exact(k))
    end
    @test length(rel) > 30
    @test maximum(rel) < 0.1
end
