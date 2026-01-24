using Test
using LinearAlgebra
using SparseArrays

using Scattensor

@testset "Matrix backend" begin
    @testset "get_length_from_localdim" begin
        @test get_length_from_localdim(16, 2) == 4
        @test get_length_from_localdim(ComplexF64[1, 0, 0, 0], 2) == 2
        @test get_length_from_localdim(Matrix(I, 8, 8), 2) == 3
    end

    @testset "operator_identity (matrix)" begin
        @test operator_identity(Matrix, 3) == Matrix{Int}(I, 3, 3)
        @test operator_identity(SparseMatrixCSC, 4) == sparse(I, 4, 4)
    end

    @testset "translation/reflection operators (matrix)" begin
        d = 2
        L = 4
        n = d^L

        T = operator_translation(SparseMatrixCSC, d, L)
        R = operator_reflection(SparseMatrixCSC, d, L)
        Iₙ = sparse(I, n, n)

        # Permutation/unitarity checks
        @test T * T' == Iₙ
        @test R * R == Iₙ
        @test T^L == Iₙ

        # Symmetry relation: T*R == R*T†
        @test T * R == R * T'
    end

    @testset "summation_local (matrix)" begin
        # Two-site local operator, compare PBC vs OBC translation invariance
        d = 2
        L = 5
        σx = SparseMatrixCSC([0 1; 1 0])
        A0 = kron(σx, σx) # acts on L0 = 2 sites

        T = operator_translation(SparseMatrixCSC, d, L)
        Hpbc = summation_local(A0, d, L; pbc = true)
        Hobc = summation_local(A0, d, L; pbc = false)

        @test norm(Hpbc * T - T * Hpbc) ≈ 0.0
        @test norm(Hobc * T - T * Hobc) > 0.0
    end

    @testset "product_locals (matrix)" begin
        d = 2
        L = 3
        σx = SparseMatrixCSC([0 1; 1 0])
        I₂ = SparseMatrixCSC(Matrix{Int}(I, 2, 2))

        # Position wraps modulo L: pos=4 -> 1
        O = product_locals(L, (σx, 4))
        @test O == kron(σx, I₂, I₂)
    end

    @testset "mathematica_format" begin
        @test mathematica_format([1 2; 3 4]) == "{{1, 2}, {3, 4}}"
    end

    @testset "fidelity_uhlmann" begin
        ρ = [0.5 0.0; 0.0 0.5]
        σ = [0.5 0.0; 0.0 0.5]
        @test isapprox(fidelity_uhlmann(ρ, σ), 1.0; atol = 1e-12)

        ρ0 = [1.0 0.0; 0.0 0.0] # |0><0|
        ρ1 = [0.0 0.0; 0.0 1.0] # |1><1|
        @test isapprox(fidelity_uhlmann(ρ0, ρ1), 0.0; atol = 1e-12)
    end
end

