"""
    periodic_even_interpolant(K, E) -> f

Build an even, 2π-periodic trigonometric interpolant from half-BZ samples.
- K: momenta in [0, π], uniform multiples of 2π/L with K[1]≈0 (π may or may not be included)
- E: energies at those K

Returns a callable f(k) that interpolates the data and satisfies f(k)=f(-k).
"""
function _periodic_even_interpolant(K, E; derivative = 0)
    @assert length(K) == length(E) && length(K) ≥ 2
    # sort by K and keep paired order
    p  = sortperm(K)
    K  = Float64.(K[p])
    E  = E[p]

    # quick sanity checks (uniform grid on [0, π], starting at 0)
    @assert K[1] ≈ 0.0 "First K must be 0"
    Δ = K[2] - K[1]
    @assert all(abs.(diff(K) .- Δ) .≤ max(1e-10, 10*eps(Float64))) "K must be uniformly spaced"
    @assert all((0.0 .- 1e-12) .≤ K .≤ (π + 1e-12)) "K must lie in [0, π]"

    M = length(K)
    # Cosine-Vandermonde: C[i,m] = cos(m * K[i]), m=0..M-1
    m = collect(0:M-1)'
    C = cos.(K .* m)
    A = C \ E  # cosine coefficients

    realE = eltype(E) <: Real

    # even, 2π-periodic interpolant
    # f_scalar(k) = (val = dot(A, cos.((0:M-1) .* k)); realE ? real(val) : val)
    # k -> isa(k, Number) ? f_scalar(k) : map(f_scalar, k)

    # precompute things we need for derivatives
    ms     = collect(0:M-1)                    # vector [0,1,2,...]
    pow_ms = derivative == 0 ? ones(eltype(A), M) : (ms .^ derivative)
    phase  = derivative * (π/2)

    # generic n-th derivative evaluator
    f_scalar = function (k)
        # cos(m k + n π/2)
        vals = cos.(ms .* k .+ phase)
        val  = dot(A, pow_ms .* vals)
        realE ? real(val) : val
    end

    k -> isa(k, Number) ? f_scalar(k) : map(f_scalar, k)
end


function dispersion_interpolant(states; derivative = 0)
    @assert !isempty(states) "No states provided"

    K = [π * Float64(s.koverpi) for s in states]
    E = [energy(s) for s in states]

    _periodic_even_interpolant(K, E, derivative = derivative)
end

export dispersion_interpolant