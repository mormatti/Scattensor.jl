# =============================================================================
# fits.jl — resonance-fitting utilities (model-independent numerics)
#
# The workhorse is the profile-linear arctan fit used throughout the
# digitization/spectroscopy projects: for FIXED nonlinear parameters
# (center x0, width w) the model
#
#     y(x) = Σ_p a_p x^p  +  A · atan( 2 (x − x0) / w )
#
# is LINEAR in (a_p, A) → exact least squares per grid node; a grid over
# (x0, w) then finds the global optimum with no degenerate local minima.
# =============================================================================

"""
    wrap_mod_pi(x) -> Float64

Map `x` to the branch `(-π/2, π/2]` by removing integer multiples of π —
the natural comparison metric for phase shifts defined mod π.
"""
wrap_mod_pi(x::Real) = x - π * round(x / π)

"""
    ArctanFit

Result of [`fit_arctan`](@ref).

# Fields
- `x0, w, A`: resonance center, full width, amplitude (total rise = `A·π`).
- `bg::Vector{Float64}`: background polynomial coefficients (order 0 first).
- `rms::Float64`: root-mean-square residual of the fit.
"""
struct ArctanFit
    x0::Float64
    w::Float64
    A::Float64
    bg::Vector{Float64}
    rms::Float64
end

"""
    fit_arctan(xs, ys; bg_order=1, x0grid=nothing, wgrid=nothing) -> ArctanFit

Fit `ys(xs)` with a polynomial background of order `bg_order` plus one arctan
resonance `A·atan(2(x−x0)/w)`, free amplitude. For each `(x0, w)` on the grid
the linear parameters are solved exactly (LSQ); the grids default to the data
range for `x0` and to logspace between `Δx/10` and the full range for `w`.

For a full Breit–Wigner resonance the rise is `A·π = π` (A = 1); a
significantly larger fitted rise signals overlapping resonances.
"""
function fit_arctan(xs::AbstractVector, ys::AbstractVector;
                    bg_order::Integer = 1,
                    x0grid = nothing, wgrid = nothing)
    n = length(xs)
    n == length(ys) || error("xs and ys must have equal length")
    n >= bg_order + 3 || error("not enough points ($n) for bg_order=$bg_order + arctan")
    lo, hi = extrema(xs)
    span = hi - lo
    x0s = x0grid === nothing ? range(lo + 0.02span, hi - 0.02span; length = 120) : x0grid
    ws  = wgrid  === nothing ? exp10.(range(log10(span / 10n), log10(span); length = 80)) : wgrid
    xm = (lo + hi) / 2
    P = [xi^p for xi in (xs .- xm), p in 0:bg_order]
    best = (Inf, Float64[], NaN, NaN)
    for x0 in x0s, w in ws
        V = hcat(P, atan.(2 .* (xs .- x0) ./ w))
        coef = V \ ys
        loss = sum(abs2, V * coef .- ys)
        loss < best[1] && (best = (loss, coef, x0, w))
    end
    loss, coef, x0, w = best
    ArctanFit(x0, w, coef[end], coef[1:end-1], sqrt(loss / n))
end
