# =============================================================================
# Luscher.jl  --  two-particle matching, phase-shift extraction, global fit
# =============================================================================
#
# Core formulas (kept visible per the assignment):
#
#     ε(k)            = E_1(k) - E_0
#     E_free(K,q)     = ε(K/2 + q) + ε(K/2 - q)
#     q L + 2 δ_K(q)  = 2π m                       (1D Lüscher quantization)
#     δ_K(q)          = π m - q L / 2
#     S_K(q)          = exp(2 i δ_K(q))
#
# K  = total momentum (= 2π Kint / L),  q = relative momentum,
# m  = Lüscher integer (free-level label).
# =============================================================================

using Optim

# ----------------------------------------------------------------------------
# Free two-particle energy:  E_free(K,q) = ε(K/2+q) + ε(K/2-q)
# ----------------------------------------------------------------------------
E_free(fit::FourierDispersion, K::Real, q::Real) =
    epsilon(fit, K / 2 + q) + epsilon(fit, K / 2 - q)

# generic bracketed root finder: all sign-change roots of f on [a,b]
function find_roots(f, a::Real, b::Real; n::Int = 2000, tol::Real = 1e-12)
    roots = Float64[]
    xs = range(a, b; length = n + 1)
    fprev = f(xs[1])
    for i in 2:length(xs)
        fcur = f(xs[i])
        if fprev == 0.0
            push!(roots, xs[i-1])
        elseif fprev * fcur < 0
            lo, hi = xs[i-1], xs[i]
            flo = fprev
            for _ in 1:80
                mid = 0.5 * (lo + hi)
                fmid = f(mid)
                if abs(fmid) < tol || (hi - lo) < tol
                    lo = hi = mid
                    break
                end
                if flo * fmid < 0
                    hi = mid
                else
                    lo, flo = mid, fmid
                end
            end
            push!(roots, 0.5 * (lo + hi))
        end
        fprev = fcur
    end
    return roots
end

# ----------------------------------------------------------------------------
# Free two-particle predictions:  q_m^(0) = π (2m − s) / L ,  s = Kint mod 2
#
# Two free particles on the ring have momenta k₁,₂ = 2π j₁,₂/L with j₁+j₂ = Kint.
# The relative momentum is q = (k₁−k₂)/2 = π(j₁−j₂)/L, and j₁−j₂ must have the
# SAME parity as Kint.  Hence for EVEN Kint the free relative momenta are
# 2πm/L, but for ODD Kint they are shifted to π(2m−1)/L (the moving-frame /
# Neveu–Schwarz half-integer shift).  s = Kint mod 2 encodes this.
# ----------------------------------------------------------------------------
"""
    free_two_particle_energies(fit, L, Kint; mrange) -> Vector{NamedTuple}

Noninteracting two-particle levels at total momentum K = 2π Kint/L: free relative
momentum q_m^(0) = π(2m − s)/L (s = Kint mod 2) and energy E_free(K, q_m^(0)).
Returns entries `(m, q0, Efree)`.
"""
function free_two_particle_energies(fit::FourierDispersion, L::Int, Kint::Int;
                                    mrange = 1:(L ÷ 2))
    K = 2π * Kint / L
    s = mod(Kint, 2)                                  # Bethe parity of the frame
    out = NamedTuple[]
    for m in mrange
        q0 = π * (2m - s) / L                         # even K: 2πm/L ; odd K: π(2m−1)/L
        (q0 < 1e-9 || q0 > π + 1e-9) && continue
        push!(out, (m = m, q0 = q0, Efree = E_free(fit, K, q0)))
    end
    return out
end

# ----------------------------------------------------------------------------
# Two-particle level matching
# ----------------------------------------------------------------------------
struct TwoParticleMatch
    L::Int
    Kint::Int
    K::Float64
    level_index::Int
    energy::Float64
    deltaE::Float64
    m::Int
    free_energy_guess::Float64
    mismatch::Float64
end

"""
    match_two_particle_levels(levels, E0, fit, L, Kint; tolerance, mmax) -> Vector{TwoParticleMatch}

Match numerical excitation energies ΔE_n = E_n - E0 in sector `Kint` to the
nearest free two-particle prediction E_free(K, q_m^(0)).  Matching is by energy
proximity and ordering.

Scientific warnings are printed when:
  * several numerical levels collide on one free level (or vice versa),
  * the mismatch exceeds `tolerance`,
  * a level lies below the two-particle threshold (possible bound state),
  * a level lies high in the spectrum (possible inelastic contamination).
"""
function match_two_particle_levels(levels::Vector{Level}, E0::Float64,
                                   fit::FourierDispersion, L::Int, Kint::Int;
                                   tolerance::Float64 = 0.15,
                                   mmax::Int = L ÷ 2,
                                   verbose::Bool = true)
    K = 2π * Kint / L
    sec = levels_in_sector(levels, Kint)                  # energy-sorted

    # two-particle continuum threshold at this K (min over relative momentum)
    qgrid = range(0, π; length = 400)
    threshold = minimum(E_free(fit, K, q) for q in qgrid)
    # 3-particle threshold ≈ 3 m1 = 1.5 × (2 m1); above it the elastic single-
    # channel extraction may be contaminated.
    inelastic = 1.5 * threshold

    # All excitations above the sector ground state (vacuum at K=0, single
    # particle at K≠0).  We classify them relative to the two-particle threshold.
    skip = 1
    excitations = [(idx = lv.index, ΔE = lv.energy - E0) for lv in sec if lv.index > skip]

    # Sub-threshold excitations are single-particle / bound states, NOT elastic
    # two-particle scattering states — exclude them and warn.  Examples: the E8
    # A2 particle (m2 = 1.618 m1 < 2 m1), or the XXZ two-magnon bound state
    # (below 2ε(0)).  A genuine K=0 scattering state has ΔE ≥ 2ε(0) = threshold,
    # equivalently E_free(K,q)=ΔE has a real relative-momentum root.
    cut = threshold * (1 - 1e-6)
    subthreshold = filter(nm -> nm.ΔE < cut, excitations)
    for nm in subthreshold
        verbose && @warn "L=$L K=$Kint level idx=$(nm.idx): ΔE=$(round(nm.ΔE,sigdigits=4)) below " *
                         "2-particle threshold $(round(threshold,sigdigits=4)); " *
                         "single-particle / BOUND STATE — excluded from elastic matching."
    end

    # Genuine two-particle candidates, energy-ordered.
    numeric = sort(filter(nm -> nm.ΔE ≥ cut, excitations), by = nm -> nm.ΔE)

    # Free predictions, energy-ordered (for K=0 this is just m = 1, 2, 3, …).
    frees = sort(free_two_particle_energies(fit, L, Kint; mrange = 1:mmax), by = fr -> fr.Efree)

    # ORDINAL matching: the n-th two-particle level ↔ the n-th free level.  This
    # is the assignment's "match by energy ordering"; proximity (mismatch) is a
    # diagnostic, not the assignment rule.  Nearest-energy matching fails for the
    # half-integer (Neveu–Schwarz) momentum quantization of free fermions.
    matches = TwoParticleMatch[]
    npairs = min(length(numeric), length(frees))
    for n in 1:npairs
        nm = numeric[n]
        fr = frees[n]
        mismatch = nm.ΔE - fr.Efree

        if abs(mismatch) > tolerance && verbose
            @warn "L=$L K=$Kint m=$(fr.m): mismatch $(round(mismatch,sigdigits=3)) > tol $tolerance " *
                  "(ΔE=$(round(nm.ΔE,sigdigits=4)) vs free $(round(fr.Efree,sigdigits=4)))."
        end
        if nm.ΔE > inelastic && verbose
            @warn "L=$L K=$Kint m=$(fr.m): ΔE=$(round(nm.ΔE,sigdigits=4)) in likely INELASTIC region " *
                  "(> ~3-particle threshold $(round(inelastic,sigdigits=4)))."
        end

        push!(matches, TwoParticleMatch(L, Kint, K, nm.idx, E0 + nm.ΔE, nm.ΔE,
                                        fr.m, fr.Efree, mismatch))
    end

    if verbose && length(numeric) != length(frees)
        @warn "L=$L K=$Kint: $(length(numeric)) two-particle candidates vs $(length(frees)) free " *
              "levels; matched $npairs (ordering may be ambiguous at the band edges)."
    end
    return matches
end

# ----------------------------------------------------------------------------
# Relative-momentum extraction:  solve ΔE = ε(K/2+q) + ε(K/2-q)
# ----------------------------------------------------------------------------
"""
    solve_relative_momentum(fit, K, deltaE; q_initial, q_bounds) -> (qbest, roots)

Solve ΔE = E_free(K,q) for the relative momentum q.  Returns the root closest to
`q_initial` together with all roots found in `q_bounds`.  Warns on multiplicity.
"""
function solve_relative_momentum(fit::FourierDispersion, K::Real, deltaE::Real;
                                 q_initial::Real = 0.0, q_bounds = (0.0, π),
                                 verbose::Bool = true)
    f(q) = E_free(fit, K, q) - deltaE
    roots = find_roots(f, q_bounds[1], q_bounds[2])
    if isempty(roots)
        verbose && @warn "No relative-momentum root for ΔE=$(round(deltaE,sigdigits=4)) at K=$(round(K,digits=3)); " *
                         "returning closest point by minimization."
        # fall back: minimize |f|
        g(q) = abs(f(q[1]))
        res = optimize(g, [clamp(q_initial, q_bounds[1], q_bounds[2])], NelderMead())
        return Optim.minimizer(res)[1], Float64[]
    end
    if length(roots) > 1 && verbose
        @warn "Multiple relative-momentum roots $(round.(roots,digits=4)) for ΔE=$(round(deltaE,sigdigits=4)); " *
              "choosing the one closest to q0=$(round(q_initial,digits=4))."
    end
    qbest = roots[argmin(abs.(roots .- q_initial))]
    return qbest, roots
end

# ----------------------------------------------------------------------------
# Phase-shift extraction:  δ_K(q) = π m − π s/2 − q L / 2 ,  s = Kint mod 2,
# S = e^{2iδ}.  (At even K, s=0, this is the textbook δ = π m − qL/2.)
# ----------------------------------------------------------------------------
struct PhaseShiftPoint
    L::Int
    Kint::Int
    K::Float64
    m::Int
    q::Float64
    delta::Float64
    S_real::Float64
    S_imag::Float64
    energy::Float64
    mismatch::Float64
end

"""
    extract_phase_shift(match, q) -> PhaseShiftPoint

δ_K(q) = π m − π s/2 − q L / 2 (s = Kint mod 2),   S = e^{2 i δ}.
The −πs/2 is the odd-frame correction; it vanishes for even Kint.
"""
function extract_phase_shift(match::TwoParticleMatch, q::Real)
    s = mod(match.Kint, 2)
    δ = π * match.m - π * s / 2 - q * match.L / 2
    S = exp(2im * δ)
    return PhaseShiftPoint(match.L, match.Kint, match.K, match.m, q, δ,
                           real(S), imag(S), match.energy, match.mismatch)
end

"""
    unwrap_phase_shifts(points) -> Vector{PhaseShiftPoint}

δ is defined only modulo π.  Group by K, sort by q, and add multiples of π to
make δ(q) continuous.  S is recomputed (it is invariant under δ → δ+π).
"""
function unwrap_phase_shifts(points::Vector{PhaseShiftPoint})
    out = PhaseShiftPoint[]
    for Kint in sort(unique(p.Kint for p in points))
        grp = sort(filter(p -> p.Kint == Kint, points), by = p -> p.q)
        isempty(grp) && continue
        δprev = grp[1].delta
        for (i, p) in enumerate(grp)
            δ = p.delta
            if i > 1
                while δ - δprev > π/2;  δ -= π;  end
                while δ - δprev < -π/2; δ += π;  end
            end
            δprev = δ
            S = exp(2im * δ)
            push!(out, PhaseShiftPoint(p.L, p.Kint, p.K, p.m, p.q, δ,
                                       real(S), imag(S), p.energy, p.mismatch))
        end
    end
    return out
end

# ----------------------------------------------------------------------------
# Wigner time delay  τ(E) = dΘ/dE = 2 dδ/dE   (S = e^{iΘ} = e^{2iδ})
#
# The two scattering particles spend an extra time τ in the interaction region.
# We give it as a function of the relative momentum q (with E the two-particle
# energy).  For free particles (δ = const) τ = 0.
# ----------------------------------------------------------------------------
"""
    time_delay_points(points, fit) -> Vector{NamedTuple}

Wigner time delay from the extracted points, by a discrete CENTRAL derivative of
δ with respect to the two-particle energy E = 2 ε(q) (K=0):
τ_i = 2 (δ_{i+1} − δ_{i−1}) / (E_{i+1} − E_{i−1}).  Points are sorted by the
relative momentum `q` (NOT the absolute energy, which is offset differently per
L), and the energy spacing is taken from the dispersion fit `ε` so that points
from different `L` combine on a single δ(E) curve.  Returns `(q, tau)`.
"""
function time_delay_points(points::Vector{PhaseShiftPoint}, fit::FourierDispersion)
    out = NamedTuple[]
    for Kint in sort(unique(p.Kint for p in points))
        grp = sort(filter(p -> p.Kint == Kint, points), by = p -> p.q)
        for i in 2:length(grp)-1
            dδ = grp[i+1].delta - grp[i-1].delta
            dE = 2 * (epsilon(fit, grp[i+1].q) - epsilon(fit, grp[i-1].q))   # ΔE = 2 ε(q)
            abs(dE) < 1e-9 && continue
            push!(out, (Kint = Kint, q = grp[i].q, tau = 2 * dδ / dE))
        end
    end
    return out
end

"""
    time_delay_from_fit(deltafun, qs, points) -> Vector{NamedTuple}

Time delay from the derivative of a smooth phase-shift fit `deltafun(q)`
(e.g. the global polynomial), evaluated at the data momenta, with dE/dq taken
from finite differences of `deltafun`'s companion energy — here we pass the
measured (q, E) and interpolate dE/dq locally.  Simpler variant: differentiate
δ(q) analytically and divide by dE/dq supplied by `dEdq`.
"""
function time_delay_from_fit(deltafun, dEdq, qs; h::Float64 = 1e-4)
    out = NamedTuple[]
    for q in qs
        dδ = (deltafun(q + h) - deltafun(q - h)) / (2h)
        push!(out, (q = q, tau = 2 * dδ / dEdq(q)))
    end
    return out
end

"""
    time_delay_exact(deltafun, epsfun, qs) -> (qs, taus)

Exact Wigner time delay τ(q) = 2 dδ/dE with E = 2 ε(q) (K=0), computed from fine
central differences of the supplied exact δ(q) and ε(q).
"""
function time_delay_exact(deltafun, epsfun, qs; h::Float64 = 1e-5)
    taus = Float64[]
    for q in qs
        dδ = (deltafun(q + h) - deltafun(q - h)) / (2h)
        dE = 2 * (epsfun(q + h) - epsfun(q - h)) / (2h)
        push!(taus, 2 * dδ / dE)
    end
    return collect(qs), taus
end

"""
    keep_lowest_per_sector(points, n) -> Vector{PhaseShiftPoint}

Keep only the `n` lowest-energy phase points in each (L, Kint) group.  The
lowest two-particle levels are the least contaminated by finite-size and
multi-particle effects, so they give the cleanest Lüscher points.
"""
function keep_lowest_per_sector(points::Vector{PhaseShiftPoint}, n::Int)
    out = PhaseShiftPoint[]
    groups = unique((p.L, p.Kint) for p in points)
    for (L, Kint) in groups
        grp = sort(filter(p -> p.L == L && p.Kint == Kint, points), by = p -> p.energy)
        append!(out, grp[1:min(n, length(grp))])
    end
    return out
end

# ----------------------------------------------------------------------------
# Convenience: full level-by-level pipeline for one (L, Kint)
# ----------------------------------------------------------------------------
"""
    phase_points_for_sector(levels, E0, fit, L, Kint; kwargs...) -> Vector{PhaseShiftPoint}

Run match → solve-q → extract-δ for one momentum sector.
"""
function phase_points_for_sector(levels::Vector{Level}, E0::Float64,
                                 fit::FourierDispersion, L::Int, Kint::Int;
                                 tolerance::Float64 = 0.15, mmax::Int = L ÷ 2,
                                 verbose::Bool = true)
    K = 2π * Kint / L
    matches = match_two_particle_levels(levels, E0, fit, L, Kint;
                                        tolerance = tolerance, mmax = mmax, verbose = verbose)
    s = mod(Kint, 2)
    pts = PhaseShiftPoint[]
    for mt in matches
        q0 = π * (2 * mt.m - s) / L
        q, roots = solve_relative_momentum(fit, K, mt.deltaE; q_initial = q0, verbose = verbose)
        if isempty(roots)        # ΔE below threshold ⇒ bound/threshold state, not scattering
            verbose && @warn "L=$L K=$Kint m=$(mt.m): no real relative-momentum root " *
                             "(ΔE below threshold); skipping as bound/threshold state."
            continue
        end
        push!(pts, extract_phase_shift(mt, q))
    end
    return pts
end

# ============================================================================
#  Global spectral fit:  δ_K(q;θ) = a0 + a1 q + a2 q^2 + a3 q^3
# ============================================================================

# polynomial phase-shift ansatz
delta_poly(θ::AbstractVector, q::Real) = sum(θ[i+1] * q^i for i in 0:length(θ)-1)

"""
    predicted_levels(fit, L, Kint, m, θ) -> (q, Epred)

For a trial phase shift θ, solve the quantization condition
q L + 2 δ_K(q;θ) = 2π m for q (near q0=2π m/L) and predict
E_pred = ε(K/2+q) + ε(K/2-q).
"""
function predicted_levels(fit::FourierDispersion, L::Int, Kint::Int, m::Int, θ::AbstractVector)
    K = 2π * Kint / L
    s = mod(Kint, 2)
    q0 = π * (2m - s) / L
    g(q) = q * L + 2 * delta_poly(θ, q) - π * (2m - s)     # qL + 2δ = π(2m − s)
    roots = find_roots(g, 1e-6, π)
    q = isempty(roots) ? q0 : roots[argmin(abs.(roots .- q0))]
    return q, E_free(fit, K, q)
end

"""
    spectral_mismatch(matches, fit, θ; weights) -> D(θ)

D(θ) = Σ_α w_α [ ΔE^num_α - E^pred_α(θ) ]^2  over matched two-particle levels.
"""
function spectral_mismatch(matches::Vector{TwoParticleMatch}, fit::FourierDispersion,
                           θ::AbstractVector; weights = nothing)
    D = 0.0
    for (i, mt) in enumerate(matches)
        _, Epred = predicted_levels(fit, mt.L, mt.Kint, mt.m, θ)
        ΔEpred = Epred - 0.0           # E_free is already an excitation energy (uses ε)
        w = weights === nothing ? 1.0 : weights[i]
        D += w * (mt.deltaE - ΔEpred)^2
    end
    return D
end

"""
    fit_phase_shift_global(matches, fit; order=3, θ0) -> (θ, cost)

Minimize the spectral mismatch D(θ) over the polynomial coefficients θ.
"""
function fit_phase_shift_global(matches::Vector{TwoParticleMatch}, fit::FourierDispersion;
                                order::Int = 3, θ0 = nothing)
    θinit = θ0 === nothing ? zeros(order + 1) : collect(θ0)
    obj(θ) = spectral_mismatch(matches, fit, θ)
    res = optimize(obj, θinit, NelderMead(), Optim.Options(iterations = 5000))
    return Optim.minimizer(res), Optim.minimum(res)
end

"""
    scan_delta_1d(matches, fit, idx, values; base) -> (values, negD)

Scan one coefficient `θ[idx]` over `values`, returning -D(θ) so that peaks mark
good spectral agreement.  `base` provides the (fixed) other coefficients.
"""
function scan_delta_1d(matches::Vector{TwoParticleMatch}, fit::FourierDispersion,
                       idx::Int, values; base = nothing)
    base = base === nothing ? zeros(maximum((idx, 2))) : collect(base)
    negD = Float64[]
    for v in values
        θ = copy(base); θ[idx] = v
        push!(negD, -spectral_mismatch(matches, fit, θ))
    end
    return collect(values), negD
end

"""
    scan_delta_2d(matches, fit, (i,j), vi, vj; base) -> (vi, vj, negD)

Two-parameter scan of θ[i], θ[j]; returns the -D(θ) grid (rows=vi, cols=vj).
"""
function scan_delta_2d(matches::Vector{TwoParticleMatch}, fit::FourierDispersion,
                       ij::Tuple{Int,Int}, vi, vj; base = nothing)
    i, j = ij
    base = base === nothing ? zeros(maximum((i, j, 2))) : collect(base)
    G = zeros(length(vi), length(vj))
    for (a, x) in enumerate(vi), (b, y) in enumerate(vj)
        θ = copy(base); θ[i] = x; θ[j] = y
        G[a, b] = -spectral_mismatch(matches, fit, θ)
    end
    return collect(vi), collect(vj), G
end
