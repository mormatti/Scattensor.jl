# Two-particle machinery on the XXZ ferromagnet — everything is exact there:
# the vacuum is a product state, σ⁻ is the exact magnon creator (ℓ = 1), the
# two-magnon sector is closed, and interactions are strictly short-range.
using SparseArrays, LinearAlgebra

@testset "pair bands" begin
    d = 2
    Jxy, Jz = 1.0, 1.5
    sp = sparse([0.0 1.0; 0.0 0.0]); sm = sparse([0.0 0.0; 1.0 0.0])
    σz = sparse([1.0 0.0; 0.0 -1.0])
    H0 = -0.25Jz * kron(σz, σz) - 0.5Jxy * (kron(sp, sm) + kron(sm, sp))
    L = 14
    H = summation_local(H0, d, L; pbc = true)
    Ω = zeros(ComplexF64, 2^L); Ω[1] = 1                 # all-up product vacuum, exact

    # --- one-body bulk functions: exact values
    S1, Hc1 = bulk_couplings(H, Ω, [Matrix(sm)], d; rmax = 2)
    s1 = S1[1, 1, :]; h1 = Hc1[1, 1, :]
    @test abs(s1[3] - 1) < 1e-12                          # s(0) = 1
    @test abs(s1[4]) < 1e-12                              # s(1) = 0 (orthogonal flips)
    @test abs(h1[3] - Jz) < 1e-12                         # h(0) = Jz (magnon mass)
    @test abs(h1[4] + Jxy / 2) < 1e-12                    # h(±1) = −Jxy/2 (hopping)

    # --- pair data
    S2, H2 = pair_couplings(H, Ω, Matrix(sm), d; seps = 1:5, hops = 4)
    # adjacent flips: one shared bond intact → measured 2Jz − Jz on the diagonal
    @test abs(H2[1, 1, 5] - Jz) < 1e-12

    # --- K = π: flat continuum at exactly 2Jz, bound state at exactly Jz
    E, isb, _ = two_particle_spectrum(π, s1, h1, S2, H2, 1:5; hops = 4, Rmax = 200)
    @test abs(minimum(E) - Jz) < 1e-8
    @test isb[argmin(E)]
    scat = E[.!isb]
    @test maximum(abs.(scat .- 2Jz)) < 1e-8               # collapsed continuum
    @test count(isb) == 1

    # --- generic K: bound level matches spin-chain ED in the same sector
    T = operator_translation(SparseMatrixCSC, d, L)
    levels = sector_spectrum(H, T, L; nlevels = 12, fullzone = true)
    E0ed = minimum(energy(s) for s in levels)
    Kint = 4                                              # K = 2π·4/14 = 4π/7
    K = 2π * Kint / L
    Em, isbm, _ = two_particle_spectrum(K, s1, h1, S2, H2, 1:5; hops = 4, Rmax = 200)
    Eb = minimum(Em[isbm])
    ed = sort([energy(s) - E0ed for s in levels
               if round(Int, abs(Float64(s.koverpi)) * L / 2) == Kint])
    @test minimum(abs.(ed .- Eb)) < 1e-2                  # bound level present in ED
    # continuum edges vs free kinematics of the measured band
    eps1(k) = Jz - Jxy * cos(k)
    lo = minimum(eps1(K / 2 + q) + eps1(K / 2 - q) for q in range(0, π; length = 400))
    hi = maximum(eps1(K / 2 + q) + eps1(K / 2 - q) for q in range(0, π; length = 400))
    scatm = Em[.!isbm]
    @test abs(minimum(scatm) - lo) < 0.02                 # wall discretization only
    @test abs(maximum(scatm) - hi) < 0.02
end
