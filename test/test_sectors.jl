# Symmetry-sector decomposition: exactness against full diagonalization.
using SparseArrays, LinearAlgebra

@testset "sectors" begin
    # Transverse-field Ising, paramagnetic phase
    d, L = 2, 8
    σx = sparse([0.0 1.0; 1.0 0.0]); σz = sparse([1.0 0.0; 0.0 -1.0])
    id2 = sparse(1.0I, 2, 2)
    H0 = -kron(σz, σz) - 2.0 * (kron(σx, id2) + kron(id2, σx)) / 2
    H = summation_local(H0, d, L; pbc = true)
    T = operator_translation(SparseMatrixCSC, d, L)

    secs = momentum_basis(T, L; fullzone = true)
    # completeness: sector dimensions sum to the full dimension
    @test sum(size(s.basis, 2) for s in secs) == d^L
    # isometry and symmetry-eigenspace property
    for s in secs
        B = s.basis
        @test norm(B' * B - I) < 1e-10
        @test norm(T * B - exp(im * π * s.label) * B) < 1e-10
    end
    # the union of all block spectra is the full spectrum (exact)
    allE = Float64[]
    for s in secs
        Hk = sector_hamiltonian(H, s)
        append!(allE, eigvals(Hermitian(Matrix(Hk))))
    end
    @test length(allE) == d^L
    @test norm(sort(allE) - eigvals(Hermitian(Matrix(H)))) < 1e-8

    # sector_spectrum agrees with dispersion_relation's lowest levels
    states = sector_spectrum(H, T, L; nlevels = 3)
    E0 = minimum(energy(s) for s in states)
    @test E0 ≈ minimum(eigvals(Hermitian(Matrix(H)))) atol = 1e-9
    # lifted states are eigenstates of H in the full space
    st = states[3]
    ψ = wavefunction(st)
    @test norm(H * ψ - energy(st) * ψ) < 1e-8

    # a generic monomial symmetry: global spin flip (order 2)
    P = kron_power(σx, L)
    flip = symmetry_blocks(P; order = 2, fullzone = true)
    @test sum(size(s.basis, 2) for s in flip) == d^L
    allE2 = Float64[]
    for s in flip
        append!(allE2, eigvals(Hermitian(Matrix(sector_hamiltonian(H, s)))))
    end
    @test norm(sort(allE2) - eigvals(Hermitian(Matrix(H)))) < 1e-8
end
