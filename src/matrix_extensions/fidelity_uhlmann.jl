# Check and write documentation for this

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