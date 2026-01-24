"""
    _periodic_even_interpolant(K, E; derivative=0) -> f

Build an even, `2π`-periodic trigonometric interpolant from samples on `[0, π]`.

This is an internal helper used by [`dispersion_interpolant`](@ref). It fits a cosine series
to the provided samples and returns a callable `f(k)` (or an array-valued map) that satisfies
`f(k) == f(-k)` by construction.

# Arguments
- `K`: Sample momenta in `[0, π]`, assumed to be uniformly spaced and starting at 0.
- `E`: Sample values (energies) at those momenta.

# Keyword Arguments
- `derivative::Integer=0`: Order of derivative to return (0 gives the interpolant itself).

# Returns
- A callable object `f` such that `f(k)` interpolates the samples (in a least-squares/linear-solve sense
  for the cosine basis) and is even and `2π`-periodic.
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


"""
    dispersion_interpolant(states; derivative=0) -> f

Construct an even, `2π`-periodic interpolant `E(k)` from a set of Bloch states.

This helper is commonly used to obtain a smooth dispersion relation from discrete samples
(typically a single band).

# Arguments
- `states`: Iterable of `BlochState` objects providing `koverpi` and `energy`.

# Keyword Arguments
- `derivative::Integer=0`: Order of derivative to return (0 gives `E(k)`).

# Returns
- A function `f(k)` that evaluates the interpolated dispersion at momentum `k` (in radians).
"""
function dispersion_interpolant(states; derivative = 0)
    @assert !isempty(states) "No states provided"

    K = [π * Float64(s.koverpi) for s in states]
    E = [energy(s) for s in states]

    _periodic_even_interpolant(K, E, derivative = derivative)
end

export dispersion_interpolant