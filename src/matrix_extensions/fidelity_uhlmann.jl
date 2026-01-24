"""
    fidelity_uhlmann(rho, sigma; root=false, tol=1e-12) -> Real

Compute the Uhlmann fidelity between two density matrices.

Given two (approximately) Hermitian, positive semidefinite matrices `rho` and `sigma`, the Uhlmann
fidelity is defined as:

`F(rho, sigma) = (tr(sqrt(sqrt(rho) * sigma * sqrt(rho))))^2`.

This implementation symmetrizes the inputs to remove small anti-Hermitian numerical noise and uses
eigendecompositions to compute matrix square roots.

# Arguments
- `rho::AbstractMatrix`: Density matrix ρ.
- `sigma::AbstractMatrix`: Density matrix σ.

# Keyword Arguments
- `root::Bool=false`: If `true`, return `tr(sqrt(sqrt(ρ) * σ * sqrt(ρ)))` (the *root fidelity*).
  If `false`, return the squared fidelity `F(ρ,σ)`.
- `tol::Real=1e-12`: Eigenvalues below `tol` are clamped to zero to improve numerical stability.

# Returns
- A real number in `[0, 1]` in the ideal PSD case (small numerical deviations may occur if inputs
  are not valid density matrices).

# Notes
- The function assumes `rho` and `sigma` are square and of the same size.
"""

function fidelity_uhlmann(rho::AbstractMatrix, sigma::AbstractMatrix; root::Bool=false, tol::Real=1e-12)
    @assert size(rho) == size(sigma) "rho and sigma must have the same size"
    n, m = size(rho)
    @assert n == m "rho and sigma must be square"

    # Symmetrize to kill tiny anti-Hermitian parts from numerics
    ρ = (rho + rho') / 2
    σ = (sigma + sigma') / 2

    # Matrix square root of a PSD Hermitian via eigen-decomposition
    # sqrt_ρ = U * sqrt.(max.(λ,0)) * U'
    evalsρ, evecsρ = eigen(Hermitian(ρ))
    λρ = max.(evalsρ, zero(eltype(evalsρ)))
    # Zero out very small negatives and tiny positives under tol
    λρ .= ifelse.(λρ .< tol, zero(eltype(λρ)), λρ)
    sqrtλρ = sqrt.(λρ)
    sqrtρ = evecsρ * Diagonal(sqrtλρ) * evecsρ'

    # Build A = sqrtρ * σ * sqrtρ, symmetrize, then take PSD square root trace
    A = sqrtρ * σ * sqrtρ
    A = (A + A') / 2  # ensure Hermitian
    evalsA = eigen(Hermitian(A)).values
    λA = max.(evalsA, zero(eltype(evalsA)))
    λA .= ifelse.(λA .< tol, zero(eltype(λA)), λA)

    tr_sqrtA = sum(sqrt.(λA))
    return root ? tr_sqrtA : tr_sqrtA^2
end