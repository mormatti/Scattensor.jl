# =============================================================================
# Dispersion.jl  --  one-particle band and ε(k) Fourier fit
# =============================================================================
#
#     ε_L(k) = E_1(L,k) - E_0(L)          (one-particle excitation energy)
#
# We fit a periodic, even cosine series
#
#     ε(k) = c_0 + Σ_{r=1}^{rmax} c_r cos(r k)
#
# which is automatically 2π-periodic and even, ε(k)=ε(-k).  This is the
# infinite-volume / interpolated dispersion used everywhere downstream.
# =============================================================================

# wrap a momentum into (-π, π]
wrap_to_pi(k) = (x = mod(k + π, 2π) - π; x ≈ -π ? π : x)

# ----------------------------------------------------------------------------
# Fourier-series dispersion object
# ----------------------------------------------------------------------------
struct FourierDispersion
    c::Vector{Float64}     # c[r+1] is the coefficient of cos(r k), r = 0..rmax
end

rmax(fit::FourierDispersion) = length(fit.c) - 1

"""
    epsilon(fit, k) -> Float64

Evaluate the fitted dispersion ε(k) = Σ_r c_r cos(r k).  Periodic in 2π and even
by construction; momentum is first mapped into (-π, π].
"""
function epsilon(fit::FourierDispersion, k::Real)
    kk = wrap_to_pi(k)
    s = 0.0
    @inbounds for r in 0:rmax(fit)
        s += fit.c[r + 1] * cos(r * kk)
    end
    return s
end

"""
    depsilon(fit, k) -> Float64

Derivative dε/dk (group velocity), used by the relative-momentum root finder.
"""
function depsilon(fit::FourierDispersion, k::Real)
    kk = wrap_to_pi(k)
    s = 0.0
    @inbounds for r in 1:rmax(fit)
        s += -fit.c[r + 1] * r * sin(r * kk)
    end
    return s
end

"""
    fit_dispersion_fourier(kvals, epsvals; rmax) -> FourierDispersion

Least-squares fit of the cosine series to the sampled band.  If there are fewer
samples than requested harmonics the order is reduced automatically.
"""
function fit_dispersion_fourier(kvals::Vector{Float64}, epsvals::Vector{Float64}; rmax::Int = 4)
    n = length(kvals)
    @assert n == length(epsvals) "kvals and epsvals length mismatch"
    R = min(rmax, n - 1)
    A = [cos(r * kvals[i]) for i in 1:n, r in 0:R]   # n × (R+1)
    c = A \ epsvals
    return FourierDispersion(collect(c))
end

# ----------------------------------------------------------------------------
# One-particle band identification (no operator overlaps; energy structure only)
# ----------------------------------------------------------------------------
"""
    one_particle_band(levels, E0, L; verbose=true) -> (kvals, epsvals)

For each momentum sector pick the one-particle candidate ε_L(k) = E_1(L,k) - E_0.

For K ≠ 0 this is simply the lowest excitation in the sector.  At K = 0 the
lowest state is the vacuum, and in a symmetry-broken / ordered phase the *second*
state can be a near-degenerate **false vacuum** rather than a particle (it shows
up as an anomalously low point breaking band continuity — see the confinement
example).  To be robust we therefore build the band from the K ≠ 0 sectors
first, extrapolate it to k = 0, and select the K = 0 excitation that is
*continuous* with that band — automatically skipping any number of false vacua.

Emits scientific warnings when the band is not clean:
  * a K=0 state looks like a second vacuum (skipped),
  * gap to the next level is small (possible crossing / not isolated),
  * ε(k) is not smooth in k.
"""
function one_particle_band(levels::Vector{Level}, E0::Float64, L::Int;
                           verbose::Bool = true, continuity_tol::Float64 = 0.25)
    # --- band from K ≠ 0 sectors (no vacuum contamination there) -----------
    kpos = Float64[]; epspos = Float64[]; gaps = Float64[]
    for Kint in momentum_sectors(levels)
        Kint ≤ 0 && continue
        sec = levels_in_sector(levels, Kint)          # energy-sorted
        isempty(sec) && continue
        E1 = sec[1].energy
        nextE = length(sec) ≥ 2 ? sec[2].energy : Inf
        push!(kpos, 2π * Kint / L); push!(epspos, E1 - E0); push!(gaps, nextE - E1)
    end

    # --- select the K = 0 one-particle state by continuity ----------------
    kvals = copy(kpos); epsvals = copy(epspos)
    if 0 in momentum_sectors(levels)
        sec0 = levels_in_sector(levels, 0)
        exc0 = [lv.energy - E0 for lv in sec0 if lv.energy > E0 + 1e-9]   # above vacuum
        if !isempty(exc0)
            chosen = exc0[1]
            if length(kpos) ≥ 2
                R = min(length(kpos) - 1, 4)
                provfit = fit_dispersion_fourier(kpos, epspos; rmax = R)
                pred0 = epsilon(provfit, 0.0)
                j = argmin(abs.(exc0 .- pred0))        # closest to band extrapolation
                chosen = exc0[j]
                bw = max(maximum(epspos) - minimum(epspos), 1e-9)
                if verbose && abs(exc0[1] - chosen) > 1e-6
                    @warn "L=$L K=0: lowest excitation ΔE=$(round(exc0[1],sigdigits=3)) looks like a " *
                          "second (near-degenerate) VACUUM, not a particle; using the " *
                          "band-continuous state ΔE=$(round(chosen,sigdigits=3)) instead."
                end
                if verbose && abs(chosen - pred0) > continuity_tol * bw
                    @warn "L=$L K=0: one-particle point ΔE=$(round(chosen,sigdigits=3)) deviates from " *
                          "band extrapolation $(round(pred0,sigdigits=3)); band may be unreliable."
                end
            end
            push!(kvals, 0.0); push!(epsvals, chosen)
        end
    end

    # order by momentum
    p = sortperm(kvals)
    kvals, epsvals = kvals[p], epsvals[p]

    if verbose && length(epsvals) ≥ 2
        bandwidth = maximum(epsvals) - minimum(epsvals)
        scale = max(bandwidth, maximum(epsvals))
        # isolation: small gap to the next level anywhere along the K>0 band
        for (kk, g) in zip(kpos, gaps)
            if isfinite(g) && g < 0.15 * scale
                @warn "L=$L one-particle band: small gap ($(round(g,sigdigits=3))) to next level " *
                      "at k=$(round(kk,digits=3)); band may not be isolated / level crossing."
            end
        end
        # smoothness: large discrete curvature relative to the bandwidth
        for i in 2:length(epsvals)-1
            curv = abs(epsvals[i+1] - 2epsvals[i] + epsvals[i-1])
            if curv > 0.5 * max(bandwidth, 1e-9)
                @warn "L=$L one-particle band: ε(k) not smooth near k=$(round(kvals[i],digits=3)) " *
                      "(curvature $(round(curv,sigdigits=3))); possible misidentification."
            end
        end
    end

    return kvals, epsvals
end
