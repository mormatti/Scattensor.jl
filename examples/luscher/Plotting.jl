# =============================================================================
# Plotting.jl  --  all figures for the Lüscher pipeline
# =============================================================================

using Plots

# ----------------------------------------------------------------------------
# (1) Energy-momentum spectrum for a given L
# ----------------------------------------------------------------------------
function plot_spectrum(levels::Vector{Level}, L::Int; E0::Float64 = NaN)
    ks = [lv.K for lv in levels]
    Es = [isnan(E0) ? lv.energy : lv.energy - E0 for lv in levels]
    ylab = isnan(E0) ? "E" : "E - E₀"
    scatter(ks ./ π, Es; xlabel = "K/π", ylabel = ylab,
            title = "Energy–momentum spectrum, L=$L",
            legend = false, markersize = 4, markerstrokewidth = 0,
            markercolor = :steelblue, grid = true)
end

# ----------------------------------------------------------------------------
# (2) One-particle band ε_L(k) and the fitted ε(k)
# ----------------------------------------------------------------------------
function plot_dispersion(kvals, epsvals, fit::FourierDispersion;
                         exact = nothing, title = "One-particle dispersion")
    ks = range(-π, π; length = 400)
    p = plot(ks ./ π, [epsilon(fit, k) for k in ks];
             label = "Fourier fit ε(k)", lw = 2, color = :crimson,
             xlabel = "k/π", ylabel = "ε(k)", title = title, legend = :top)
    scatter!(p, kvals ./ π, epsvals; label = "ε_L(k) data",
             markersize = 5, markerstrokewidth = 0, color = :black)
    if exact !== nothing
        plot!(p, ks ./ π, [exact(k) for k in ks];
              label = "exact (free fermion)", lw = 2, ls = :dash, color = :seagreen)
    end
    return p
end

# ----------------------------------------------------------------------------
# (3) Two-particle: numerical excitations vs free predictions at fixed K
# ----------------------------------------------------------------------------
# Numerical (matched) vs free prediction at the SAME Lüscher integer m — the
# vertical gap is the interaction shift that the phase shift encodes.  Levels the
# matcher excluded (bound states / below the two-particle threshold) are shown
# separately on the left, so the m-axis is not offset by them.
function plot_two_particle_levels(numeric_ΔE::Vector{Float64},
                                  frees, matches::Vector{TwoParticleMatch}, Kint::Int, L::Int)
    p = plot(; xlabel = "Lüscher integer m", ylabel = "ΔE = E - E₀",
             title = "Two-particle levels, L=$L, Kint=$Kint", legend = :topleft)
    if !isempty(matches)
        ms  = [m.m for m in matches]
        num = [m.deltaE for m in matches]
        fre = [m.free_energy_guess for m in matches]
        for (m, a, b) in zip(ms, num, fre)             # connector = mismatch
            plot!(p, [m, m], [a, b]; color = :gray, lw = 1, label = "")
        end
        scatter!(p, ms, fre; label = "free prediction E_free(m)", color = :seagreen,
                 markershape = :diamond, markersize = 7, markerstrokewidth = 0)
        scatter!(p, ms, num; label = "numerical (matched)", color = :crimson,
                 markersize = 6, markerstrokewidth = 0)
        # excluded sub-threshold (bound) levels: numerical levels below the lowest match
        matchedE = num
        lowmatch = minimum(matchedE)
        bound = [e for e in numeric_ΔE if e < lowmatch - 1e-6 &&
                 all(abs(e - me) > 1e-6 for me in matchedE)]
        if !isempty(bound)
            scatter!(p, fill(0.5, length(bound)), bound; label = "excluded (bound / < threshold)",
                     color = :gray, markershape = :xcross, markersize = 6, markerstrokewidth = 1)
        end
    end
    return p
end

# ----------------------------------------------------------------------------
# (4) Extracted phase shift δ_K(q) vs q  (the headline plot)
# ----------------------------------------------------------------------------
# fold δ (defined mod π) to the representative nearest a reference value
_fold_pi(δ, ref) = δ + π * round((ref - δ) / π)

function plot_phase_shift(points::Vector{PhaseShiftPoint};
                          benchmark = nothing, fold_ref = nothing,
                          title = "Phase shift δ_K(q)")
    p = plot(; xlabel = "q", ylabel = "δ_K(q)", title = title, legend = :best)
    for Kint in sort(unique(pt.Kint for pt in points))
        grp = sort(filter(pt -> pt.Kint == Kint, points), by = pt -> pt.q)
        isempty(grp) && continue
        δvals = fold_ref === nothing ? [pt.delta for pt in grp] :
                                       [_fold_pi(pt.delta, fold_ref(pt.q)) for pt in grp]
        scatter!(p, [pt.q for pt in grp], δvals;
                 label = "L-data, Kint=$Kint", markersize = 6, markerstrokewidth = 0)
    end
    if benchmark !== nothing
        bq, bδ, blabel = benchmark
        plot!(p, bq, bδ; label = blabel, lw = 2, ls = :dash, color = :black)
    end
    return p
end

# ----------------------------------------------------------------------------
# (4b) Phase shift δ_K(q) at SEVERAL total momenta K — the multi-K view
# ----------------------------------------------------------------------------
# One colour per Kint sector: extracted points + (if given) the K-dependent exact
# curve δ_exact(Kint, q).  Points are folded (mod π) to their own K's curve.  At
# K=0 the curve is L-independent; at K≠0 it depends on K=2π Kint/L, so this plot
# is made at a single L (the exact curve must match that L).
function plot_phase_shift_multiK(points::Vector{PhaseShiftPoint};
                                 delta_exact = nothing, fold_ref = nothing,
                                 title = "Phase shift δ_K(q), all K")
    p = plot(; xlabel = "q", ylabel = "δ_K(q)", title = title, legend = :best)
    Kints = sort(unique(pt.Kint for pt in points))
    qcurve = range(0.03, π - 0.02; length = 250)
    for (ci, Kint) in enumerate(Kints)
        grp = sort(filter(pt -> pt.Kint == Kint, points), by = pt -> pt.q)
        isempty(grp) && continue
        col = ci
        # fold (mod π) to the exact curve if given, else to a scalar reference,
        # else show the raw unwrapped branch
        δvals = delta_exact !== nothing ? [_fold_pi(pt.delta, delta_exact(Kint, pt.q)) for pt in grp] :
                fold_ref    !== nothing ? [_fold_pi(pt.delta, fold_ref) for pt in grp] :
                                          [pt.delta for pt in grp]
        scatter!(p, [pt.q for pt in grp], δvals; label = "data, Kint=$Kint",
                 color = col, markersize = 6, markerstrokewidth = 0)
        if delta_exact !== nothing
            plot!(p, qcurve, [delta_exact(Kint, q) for q in qcurve];
                  label = "exact, Kint=$Kint", color = col, lw = 2, ls = :dash)
        end
    end
    # when there is no exact curve, mark the free reference the points are folded to
    if delta_exact === nothing && fold_ref !== nothing
        hline!(p, [fold_ref]; label = "free (δ=$(round(fold_ref, digits=3)))",
               color = :black, lw = 1.5, ls = :dot)
    end
    return p
end

# ----------------------------------------------------------------------------
# (5) Extracted 2 δ_K(q) vs q
# ----------------------------------------------------------------------------
function plot_two_delta(points::Vector{PhaseShiftPoint}; title = "2 δ_K(q)")
    p = plot(; xlabel = "q", ylabel = "2 δ_K(q)", title = title, legend = :best)
    for Kint in sort(unique(pt.Kint for pt in points))
        grp = sort(filter(pt -> pt.Kint == Kint, points), by = pt -> pt.q)
        isempty(grp) && continue
        scatter!(p, [pt.q for pt in grp], [2 * pt.delta for pt in grp];
                 label = "Kint=$Kint", markersize = 6, markerstrokewidth = 0)
    end
    return p
end

# ----------------------------------------------------------------------------
# (6) S_K(q) real and imaginary parts vs q
# ----------------------------------------------------------------------------
function plot_smatrix(points::Vector{PhaseShiftPoint}; benchmark = nothing,
                     title = "S_K(q) = e^{2iδ}")
    grp = sort(points, by = pt -> pt.q)
    qs = [pt.q for pt in grp]
    p = plot(; xlabel = "q", ylabel = "S", title = title, legend = :best)
    if benchmark !== nothing
        bq, bre, bim, blabel = benchmark
        plot!(p, bq, bre; label = "Re S ($blabel)", color = :crimson, lw = 2, ls = :dash)
        plot!(p, bq, bim; label = "Im S ($blabel)", color = :steelblue, lw = 2, ls = :dash)
    end
    scatter!(p, qs, [pt.S_real for pt in grp]; label = "Re S (extracted)", color = :crimson,
             markersize = 6, markerstrokewidth = 0)
    scatter!(p, qs, [pt.S_imag for pt in grp]; label = "Im S (extracted)", color = :steelblue,
             markersize = 6, markerstrokewidth = 0)
    return p
end

# ----------------------------------------------------------------------------
# (7) Global-fit comparison: numerical vs Lüscher-predicted levels
# ----------------------------------------------------------------------------
function plot_global_fit(matches::Vector{TwoParticleMatch}, fit::FourierDispersion,
                         θ::AbstractVector; title = "Global fit: numerical vs predicted")
    num = [m.deltaE for m in matches]
    pred = [predicted_levels(fit, m.L, m.Kint, m.m, θ)[2] for m in matches]
    p = scatter(num, pred; xlabel = "ΔE numerical", ylabel = "ΔE predicted (θ)",
                title = title, legend = :topleft, label = "matched levels",
                markersize = 5, markerstrokewidth = 0, color = :purple)
    lo, hi = minimum(num), maximum(num)
    plot!(p, [lo, hi], [lo, hi]; label = "y = x", color = :black, ls = :dash)
    return p
end

# ----------------------------------------------------------------------------
# (8) Phase-shift ansatz scans (peaks of -D mark good agreement)
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Wigner time delay τ(q) = 2 dδ/dE : exact curve + data (central differences)
# ----------------------------------------------------------------------------
function plot_time_delay(tau_points, exact_curve;
                        fit_points = nothing, title = "Wigner time delay τ(q)")
    p = plot(; xlabel = "q", ylabel = "τ(q) = 2 dδ/dE", title = title, legend = :best)
    if exact_curve !== nothing
        eq, eτ = exact_curve
        plot!(p, eq, eτ; label = "exact 2 dδ/dE", lw = 2, ls = :dash, color = :black)
    end
    if fit_points !== nothing
        plot!(p, [t.q for t in fit_points], [t.tau for t in fit_points];
              label = "from δ-fit derivative", lw = 2, color = :seagreen)
    end
    if !isempty(tau_points)
        scatter!(p, [t.q for t in tau_points], [t.tau for t in tau_points];
                 label = "from data (central diff)", markersize = 6,
                 markerstrokewidth = 0, color = :crimson)
    end
    return p
end

function plot_scan_1d(values, negD; label = "θ", title = "Spectral scan −D(θ)")
    plot(values, negD; xlabel = label, ylabel = "−D(θ)", title = title,
         lw = 2, legend = false, color = :darkorange)
end

function plot_scan_2d(vi, vj, G; xlabel = "θ_i", ylabel = "θ_j",
                      title = "Spectral scan −D(θ)")
    heatmap(vj, vi, G; xlabel = ylabel, ylabel = xlabel, title = title, color = :viridis)
end
