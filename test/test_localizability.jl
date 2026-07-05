# Localizability: exact analytic checks on hand-built Bloch states.
using LinearAlgebra

@testset "localizability" begin
    d, L = 2, 8
    up = [1.0 + 0im, 0.0]; dn = [0.0, 1.0 + 0im]
    basis(j) = j == 0 ? up : dn
    ket(bits) = kron([basis(b) for b in bits]...)
    Ω = ket(zeros(Int, L))

    # --- one-magnon Bloch state (W state at k = 0): perfect 1-site particle
    W = sum(ket([i == j ? 1 : 0 for i in 1:L]) for j in 1:L) / sqrt(L)
    for ℓ in 1:L-1
        @test localizability(W, Ω, ℓ) ≈ 1.0 atol = 1e-10
    end
    @test support_size(W, Ω) ≈ 1.0 atol = 1e-9

    # --- two-site bound cluster (adjacent flips, k = 0): C̃(1) = 0, size 2
    C2 = sum(ket([i == j || i == mod1(j + 1, L) ? 1 : 0 for i in 1:L]) for j in 1:L) / sqrt(L)
    @test transition_trace_norm(C2, Ω, 1) < 1e-12          # one site cannot create it
    for ℓ in 2:L-1
        # ℓ contiguous sites host ℓ−1 full placements of the 2-site cluster:
        # ‖A‖₁ = (ℓ−1)/√(ℓ L) exactly, hence C̃ = (ℓ−1)/ℓ … up to the wrap term
        @test localizability(C2, Ω, ℓ) >= (ℓ - 1) / ℓ - 1e-9
    end
    @test 1.5 < support_size(C2, Ω) < 3.0

    # --- generic random state: C̃ small at small ℓ (no compact creator)
    rnd = normalize!(randn(ComplexF64, d^L))
    @test localizability(rnd, Ω, 1) < 0.7
    @test support_size(rnd, Ω) > support_size(W, Ω)

    # --- optimal creator attains the trace norm and is unitary
    φ, r = optimal_creator(C2, Ω, 3)
    @test norm(φ' * φ - I) < 1e-9
    Mψ = reshape(C2, d^3, d^(L - 3)); MΩ = reshape(Ω, d^3, d^(L - 3))
    A = Matrix(Mψ * MΩ')
    @test abs(tr(A' * φ)) ≈ sum(svdvals(A)) atol = 1e-10
end
