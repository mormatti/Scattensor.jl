using KrylovKit
using ITensors, ITensorMPS

# Assumptions
# - The Hamiltonian is an Hermitian operator;
# - The system is finite, PBC and translationally invariant (by now);
# - For every value of the momentum k expressed as a phase (e^ik), k / 2π is rational
#   this assumption is quite general, and includes SSB case.

function excitations(H::Hermitian, T::Matrix, L::Int, k::Real, n::Int; kwargs...)
    #TODO: Implement the function
end

function excitations_exact(H::Hermitian, T::Matrix, L::Int, k::Real, n::Int)
    # Assertions
    N = size(H)[1]
    @assert size(T) == (N, N) "T must be a square matrix, but size(T) = $(size(T))"
    @assert H * T ≈ T * H "The Hamiltonian must be translationally invariant."
    @assert T * T' ≈ T' * T ≈ I "The translation operator must be unitary."

    # We take the modulus 2π of the momentum.
    k = k % 2π

    # We define the unit of phases.
    Tk = e^(im * k) * T

    # We project the Hamiltonian to the eigenspace of momentum k.
    Hk = H
    for _ in 1:(L-1)
        Hk += Tk * Hk
    end
    Hk /= L
    Hk = Hermitian(Hk)

    # Compute the eigenvalues of Hk and we return the first n values and vectors.
    λ, v = eigen(Hk)
    return λ[1:n], [v[:, i] for i in 1:n]
end

function _excitations_krylov(k::Real, H::Function, L::Int, T::Function, n::Int)
    # In case H and T are matrices, assert Hermiticity, unitarity
    # and translational invariance.

    # We take the modulus 2π of the momentum.
    k = k % 2π

    # We define the projected hamiltonian function.
    function Hk(v)
        v = H(v)
        for _ in 1:(L-1)
            v = e^(im * k) * T(v) + v
        end
        return v / L
    end

    # Compute the eigenvalues of Hk, and we sleect the first n values and vectors.
    energies, states, info = eigsolve(Hk, n; howmany = n, which = :SR, T = ComplexF64)

    # We return the energies and the first n states.
    return energies, states
end

# Create a funtion f(H,T,L,n) which search in the 0-momentum eigenspace of H.

# N = size(H)[1]
# @assert size(T) == (N, N) "T must be a square matrix, but size(T) = $(size(T))"
# @assert H * T ≈ T * H "The Hamiltonian must be translationally invariant."
# @assert T * T' ≈ T' * T ≈ I "The translation operator must be unitary."
function _excitations_zero_k_krylov(H::AbstractMatrix, T::AbstractMatrix, L::Int, n::Int)
    # We project the Hamiltonian to the eigenspace of momentum 0.
    H0 = H
    for _ in 1:(L-1)
        H0 += T * H0
    end
    H0 /= L

    # Compute the eigenvalues of Hk and we return the first n values and vectors.
    λ, v = eigen(Hk)
    return λ[1:n], [v[:, i] for i in 1:n]
end