# Bulk-operator band reconstruction (OBC, no translation operator needed)
# against the exact PBC dispersion — TFI deep paramagnet.
using SparseArrays, LinearAlgebra

@testset "bulk bands" begin
    d = 2
    σx = sparse([0.0 1.0; 1.0 0.0]); σz = sparse([1.0 0.0; 0.0 -1.0])
    id2 = sparse(1.0I, 2, 2)
    J, h = 1.0, 4.0
    H0 = -J * kron(σz, σz) - 0.5h * (kron(σx, id2) + kron(id2, σx))

    # --- exact reference band from PBC sector ED at L = 12
    Lp = 12
    Hp = summation_local(H0, d, Lp; pbc = true)
    Tp = operator_translation(SparseMatrixCSC, d, Lp)
    exact = sector_spectrum(Hp, Tp, Lp; nlevels = 2, fullzone = true)
    E0p = minimum(energy(s) for s in exact)
    ref = Dict{Float64, Float64}()
    for s in exact
        dE = energy(s) - E0p
        dE < 1e-9 && continue
        k = Float64(s.koverpi) * π
        ref[k] = min(get(ref, k, Inf), dE)
    end

    # --- bulk reconstruction on an OBC chain, dressed creator from the trap
    Lo = 12
    Ho = summation_local(H0, d, Lo; pbc = false)
    Fo = eigen(Hermitian(Matrix(Ho)))
    Ω = normalize!(Vector{ComplexF64}(Fo.vectors[:, 1]))
    To = operator_translation(SparseMatrixCSC, d, Lo)          # only to build the trap
    Htrap = deformed_hamiltonian(H0, d, Lo; T = To)
    _, seeds = trap_seeds(Htrap, 2)
    # dressed 3-site flip, extracted WHERE THE SEED LIVES (window scan)
    seed = Vector{ComplexF64}(normalize(seeds[2]))
    seed -= dot(Ω, seed) * Ω; normalize!(seed)
    tns = [transition_trace_norm(seed, Ω, 3; site = s) for s in 1:Lo-2]
    φ, _ = optimal_creator(seed, Ω, 3; site = argmax(tns))

    S, Hc = bulk_couplings(Ho, Ω, [φ], d; rmax = 3)
    @test abs(S[1, 1, 4]) > 0.5                                # normalized-ish at r = 0
    @test abs(S[1, 1, 1]) < abs(S[1, 1, 3])                    # overlaps decay with |r|
    ks, bands = wannier_bands(S, Hc; nk = 100)

    # compare at the exact momenta: RELATIVE error (band values are O(h))
    errs = Float64[]
    for (k, E) in ref
        i = argmin(abs.(ks .- k))
        isnan(bands[i, 1]) && continue
        push!(errs, abs(bands[i, 1] - E) / E)
    end
    @test length(errs) >= 5
    @test maximum(errs) < 0.03          # single-species frame: ≲2% across the zone
    # variance certificate: the reconstructed band is far below the 2-particle
    # continuum, so the corresponding trapped seed has small energy variance
    @test energy_variance(Ho, normalize(Vector{ComplexF64}(seeds[1]))) <
          energy_variance(Ho, normalize!(randn(ComplexF64, size(Ho, 1))))
end
