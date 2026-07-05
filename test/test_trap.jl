# Trap (SSD-deformed Hamiltonian) → effective bands vs exact sector spectrum.
using SparseArrays, LinearAlgebra

@testset "trap" begin
    d, L = 2, 8
    σx = sparse([0.0 1.0; 1.0 0.0]); σz = sparse([1.0 0.0; 0.0 -1.0])
    id2 = sparse(1.0I, 2, 2)
    J, h = 1.0, 4.0                       # deep paramagnet: isolated 1-particle band
    H0 = -J * kron(σz, σz) - 0.5h * (kron(σx, id2) + kron(id2, σx))
    H = summation_local(H0, d, L; pbc = true)
    T = operator_translation(SparseMatrixCSC, d, L)

    Hp = deformed_hamiltonian(H0, d, L; T = T)
    @test norm(Hp - Hp') < 1e-10
    @test norm(Hp * T - T * Hp) > 1e-3            # the trap breaks translations

    en, seeds = trap_seeds(Hp, 8)
    @test issorted(en)

    H̃, T̃, W = effective_pair(H, T, seeds, L)
    @test norm(H̃ - H̃') < 1e-8
    @test norm(T̃' * T̃ - I) < 1e-6                 # effective translation stays unitary

    # exact reference from the sector decomposition
    exact = sector_spectrum(H, T, L; nlevels = 3, fullzone = true)
    eff = dispersion_relation(Matrix(H̃), Matrix(T̃), L; nlevels = 3)
    E0ex = minimum(energy(s) for s in exact)
    E0eff = minimum(energy(s) for s in eff)
    @test E0eff ≈ E0ex atol = 1e-5                # variational vacuum, tight

    # one-particle band: LOWEST excitation per momentum sector, both sides
    function lowest_band(states, E0)
        m = Dict{Rational, Float64}()
        for s in states
            dE = energy(s) - E0
            dE < 1e-9 && continue
            m[s.koverpi] = min(get(m, s.koverpi, Inf), dE)
        end
        m
    end
    bex = lowest_band(exact, E0ex)
    bef = lowest_band(eff, E0eff)
    matched = [k for k in keys(bef) if haskey(bex, k) && abs(bex[k] - bef[k]) < 5e-3]
    @test length(matched) >= 4                    # band recovered across the zone
end
