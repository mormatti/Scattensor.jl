# Composite channels + ab-initio Friedrichs–Lee in the two-particle sector.
#
# Three independent checks:
#  (i)   ANALYTIC surface-impurity toy for `resonance_parameters`: a single
#        impurity level coupled to the surface site of an open tight-binding
#        chain has an exactly known width Γ = 2π g₀² ρ₁(E) = g₀² √(4−E²).
#  (ii)  XXZ end-to-end through `friedrichs_lee`: the 2-site adjacent-pair
#        creator is (proportional to) the exact bound-pair creator, so it must
#        be dropped as a Löwdin-null channel while the enlarged spectrum stays
#        identical to the pair-only spectrum (no double counting).
#  (iii) INDEX correctness: one cross table sX rebuilt by explicit dense kron.
using SparseArrays, LinearAlgebra

@testset "composite channels" begin

    # =====================================================================
    # (i) analytic surface-impurity Friedrichs–Lee toy
    # =====================================================================
    @testset "surface-impurity resonance width" begin
        N = 600                                   # open chain, hopping t = 1
        g0 = 0.3
        Echi = 0.5                                # bare level inside the band
        # open tight-binding chain: continuum [-2, 2], surface DOS
        # ρ₁(E) = √(4 − E²) / (2π); eigenpairs give εₐ and the amplitude on site 1.
        Hchain = SymTridiagonal(zeros(N), ones(N - 1))
        F = eigen(Matrix(Hchain))
        eps = F.values                            # N continuum levels in (-2, 2)
        g = g0 .* F.vectors[1, :]                 # gₐ = g₀ ⟨site 1 | a⟩

        # η tuning: 5 mean level spacings of the continuum window (the
        # `resonance_parameters` default). A few spacings are needed so the
        # Gaussian-kernel DOS estimate is smooth; the residual O(η²) curvature
        # bias is tiny here (band is smooth at E_R ≈ 0.5). Measured rel. error
        # ≈ 2e-4 — orders of magnitude below the 10% budget; η from ~1 spacing
        # up to ~0.1 (≈15 spacings) all stay < 1%.
        η = 5 * (maximum(eps) - minimum(eps)) / (N - 1)
        rp = resonance_parameters(Echi, g, eps; eta = η)

        # analytic width AT THE EXTRACTED resonance position
        Γ_exact = g0^2 * sqrt(4 - rp.E_R^2)
        relerr = abs(rp.Gamma - Γ_exact) / Γ_exact
        @test relerr < 0.10                       # (holds with a huge margin)
        @test minimum(eps) < rp.E_R < maximum(eps)  # resonance sits in the band

        # bound (below-continuum) channel: Echi = 3.0 lies outside [-2, 2]
        rb = resonance_parameters(3.0, g, eps; eta = η)
        @test rb.Gamma == 0.0                     # pole off the cut → no width
        @test rb.E_R > maximum(eps)               # stays outside the continuum
    end

    # =====================================================================
    # XXZ testbed shared by (ii) and (iii): all-up ferromagnet vacuum,
    # σ⁻ the exact ℓ = 1 magnon creator, two-magnon sector closed.
    # =====================================================================
    d = 2
    Jxy, Jz = 1.0, 1.5
    sp = sparse([0.0 1.0; 0.0 0.0]); sm = sparse([0.0 0.0; 1.0 0.0])
    σz = sparse([1.0 0.0; 0.0 -1.0])
    H0 = -0.25Jz * kron(σz, σz) - 0.5Jxy * (kron(sp, sm) + kron(sm, sp))

    # =====================================================================
    # (ii) end-to-end Friedrichs–Lee: composite = 2-site adjacent-pair creator
    # =====================================================================
    @testset "XXZ Friedrichs–Lee: drop + no double counting" begin
        L = 14
        H = summation_local(H0, d, L; pbc = true)
        Ω = zeros(ComplexF64, 2^L); Ω[1] = 1

        S1, Hc1 = bulk_couplings(H, Ω, [Matrix(sm)], d; rmax = 2)
        s1 = S1[1, 1, :]; h1 = Hc1[1, 1, :]
        S2, H2 = pair_couplings(H, Ω, Matrix(sm), d; seps = 1:5, hops = 4)

        # χ̂ = σ⁻ ⊗ σ⁻ creates the adjacent pair at (c, c+1): a genuine 2-site
        # composite creator. At any K it equals e^{-iK/2}|K, r=1⟩ (a pair-basis
        # vector, since r = 1 ∈ seps), i.e. exactly linearly dependent.
        χ = kron(Matrix(sm), Matrix(sm))
        S0, H0c = bulk_couplings(H, Ω, [χ], d; rmax = 4)
        sX, hX = composite_pair_cross(H, Ω, [χ], Matrix(sm), d; seps = 1:5, hops = 4)

        # --- K = π: composite must be DROPPED, spectrum unchanged (E_b = Jz)
        flπ = friedrichs_lee(π, S0, H0c, sX, hX, s1, h1, S2, H2, 1:5; hops = 4, Rmax = 200)
        Eπ, isbπ, _ = two_particle_spectrum(π, s1, h1, S2, H2, 1:5; hops = 4, Rmax = 200)
        @test !isempty(flπ.dropped)                        # (a) near-null channel dropped
        @test flπ.dropped == [1] && isempty(flπ.kept)
        @test minimum(abs.(flπ.energies .- Jz)) < 1e-6     # (b) bound level E_b(π) = Jz kept
        @test length(flπ.energies) == length(Eπ)           # no extra state → no double counting
        @test maximum(abs.(sort(flπ.energies) .- sort(Eπ))) < 1e-6   # continuum unchanged
        @test count(flπ.isbound) == 1

        # --- K = π/2 (non-degenerate): enlarged bound energy = pair-only bound
        K = π / 2
        flh = friedrichs_lee(K, S0, H0c, sX, hX, s1, h1, S2, H2, 1:5; hops = 4, Rmax = 200)
        Eh, isbh, _ = two_particle_spectrum(K, s1, h1, S2, H2, 1:5; hops = 4, Rmax = 200)
        @test !isempty(flh.dropped)
        Eb_fl = minimum(flh.energies[flh.isbound])
        Eb_pair = minimum(Eh[isbh])
        @test abs(Eb_fl - Eb_pair) < 1e-6                  # pair basis already complete there
        @test length(flh.energies) == length(Eh)
    end

    # =====================================================================
    # (iii) index correctness: one sX table vs explicit dense kron
    # =====================================================================
    @testset "composite_pair_cross indices" begin
        L = 8
        H = summation_local(H0, d, L; pbc = true)
        Ω = zeros(ComplexF64, 2^L); Ω[1] = 1
        χ = kron(Matrix(sm), Matrix(sm))
        hops = 1; seps = [1, 2]; c = 3
        sX, _ = composite_pair_cross(H, Ω, [χ], Matrix(sm), d; seps = seps, hops = hops, center = c)

        # independent dense embedding 1_{j-1} ⊗ op ⊗ 1_{L-j-ℓ+1}
        emb(op, j) = begin
            ℓ = round(Int, log(d, size(op, 1)))
            kron(Matrix{ComplexF64}(I, d^(j - 1), d^(j - 1)),
                 kron(ComplexF64.(op), Matrix{ComplexF64}(I, d^(L - j - ℓ + 1), d^(L - j - ℓ + 1))))
        end
        Pv(v) = v - dot(Ω, v) * Ω
        sX_man = zeros(ComplexF64, 1, length(seps), 2hops + 1)
        for (rj, r) in enumerate(seps), m in -hops:hops
            ket = Pv(emb(Matrix(sm), c + m) * (emb(Matrix(sm), c + m + r) * Ω))  # pair at (c+m, c+m+r)
            bra = Pv(emb(χ, c) * Ω)                                             # composite, left edge c
            sX_man[1, rj, m + hops + 1] = dot(bra, ket)
        end
        @test maximum(abs.(sX .- sX_man)) < 1e-10
        # the single nonzero overlap: composite {c,c+1} = pair (c, c+1) ⇒ m=0, r=1
        @test abs(sX[1, 1, hops + 1] - 1) < 1e-10
    end

end
