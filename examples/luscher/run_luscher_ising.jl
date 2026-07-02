# =============================================================================
# run_luscher_ising.jl  --  1D Lüscher phase-shift extraction, end to end
# =============================================================================
#
# Pipeline:
#   spectrum  →  one-particle band ε(k)  →  Fourier fit  →  two-particle level
#   matching  →  relative momentum q  →  phase shift δ_K(q), S=e^{2iδ}
#   →  global spectral fit  →  plots + saved results.
#
# Scenarios:
#   (A) XXZ FERROMAGNET  Δ=1.5, h=2  — HEADLINE analytic benchmark.  Bethe-ansatz
#                                      integrable: extracted δ(q)/S(q) compared to
#                                      the EXACT two-magnon S(q)=−(1−Δe^{−iq})/(1−Δe^{iq}).
#   (B) ISING FREE FERMION h_x=4, h_z=0 — exact dispersion ε(k)=2√(J²+h_x²−2Jh_x cosk);
#                                      flat phase S=−1 (δ=π/2).  Deep in the disordered
#                                      phase so the band is separated from multi-particle states.
#   (B2) ISING DISORDERED h_x=4, h_z=1 — interacting & NON-integrable, but with a
#                                      cleanly SEPARATED band (unlike critical E8): the
#                                      free-fermion point with the longitudinal field ON,
#                                      so δ_K(q) deviates smoothly from π/2.
#   (C) ISING E8 POINT   h_x=1, h_z=0.04 — literature benchmark: golden-ratio mass
#                                      ratio m₂/m₁→1.6180 (Zamolodchikov; Coldea 2010).
#   (C2) ISING E8 ISOLATED h_x=1, h_z=0.5 — same critical point, larger field: the band
#                                      isolates (clean phase) but drifts off exact E8
#                                      scaling (golden ratio 0.09%→~1.4%): a crossover.
#   (D) ISING CONFINEMENT h_x=0.6, h_z=0.05 — assignment defaults; NON-integrable
#                                      (no exact δ); demonstrates false-vacuum handling.
#
# Most scenarios run total momentum K=0 only; XXZ, free-fermion and the disordered
# point also extract nonzero-K sectors (`multiK`) → `phaseshift_multiK.png`, with
# the K-dependent exact δ_K(q) overlaid where a closed form exists.
#
# Run from the package root:
#     julia --project=. examples/luscher/run_luscher_ising.jl
# =============================================================================

using Scattensor
using SparseArrays, LinearAlgebra
using Plots
gr()

const HERE = @__DIR__
include(joinpath(HERE, "Model.jl"))
include(joinpath(HERE, "Spectrum.jl"))
include(joinpath(HERE, "Dispersion.jl"))
include(joinpath(HERE, "Luscher.jl"))
include(joinpath(HERE, "Benchmark.jl"))
include(joinpath(HERE, "Plotting.jl"))

const PLOTS   = joinpath(HERE, "plots")
const RESULTS = joinpath(HERE, "results")
mkpath(PLOTS); mkpath(RESULTS)

# ----------------------------------------------------------------------------
# minimal dependency-free JSON writer (numbers / strings / flat vectors / dicts)
# ----------------------------------------------------------------------------
jsonval(x::Real)   = isfinite(x) ? string(x) : "null"
jsonval(x::String) = "\"" * x * "\""
jsonval(x::AbstractVector) = "[" * join(jsonval.(x), ", ") * "]"
function write_json(path, d::AbstractDict)
    open(path, "w") do io
        println(io, "{")
        keys_ = collect(keys(d))
        for (i, k) in enumerate(keys_)
            comma = i < length(keys_) ? "," : ""
            println(io, "  \"$k\": ", jsonval(d[k]), comma)
        end
        println(io, "}")
    end
end

function write_phase_csv(path, pts::Vector{PhaseShiftPoint})
    open(path, "w") do io
        println(io, "L,Kint,K,m,q,delta,two_delta,S_real,S_imag,energy,mismatch")
        for p in pts
            println(io, join((p.L, p.Kint, p.K, p.m, p.q, p.delta, 2p.delta,
                              p.S_real, p.S_imag, p.energy, p.mismatch), ","))
        end
    end
end

# ============================================================================
#  One full analysis for a given parameter set
# ============================================================================
function analyze(name::String, params, config::LuscherConfig;
                 benchmark::Symbol = :none, keep_lowest::Int = 2,
                 multiK::Vector{Int} = Int[], delta_exact_K = nothing)
    println("\n" * "="^70)
    println("SCENARIO $name :  ", params)
    println("="^70)

    # per-scenario output folders (keeps plots/ and results/ tidy)
    pdir = joinpath(PLOTS, name);   mkpath(pdir)
    rdir = joinpath(RESULTS, name); mkpath(rdir)

    # --- (1) spectra for every L ------------------------------------------
    levels_by_L = Dict{Int, Vector{Level}}()
    E0_by_L     = Dict{Int, Float64}()
    for L in config.Ls
        println("\n[L=$L] building and diagonalizing ...")
        H = build_hamiltonian(L, params)
        T = build_translation(L, params)
        levels = compute_momentum_spectrum(H, T, L; nlevels = config.nlevels,
                                           parity_op = reflection_operator(L, params))
        E0, _ = vacuum_energy(levels)
        levels_by_L[L] = levels
        E0_by_L[L] = E0
        savefig(plot_spectrum(levels, L; E0 = E0),
                joinpath(pdir, "spectrum_L$(L).png"))
    end

    # --- (2) one-particle band + Fourier fit (from the largest L) ---------
    Lbig = maximum(config.Ls)
    kvals, epsvals = one_particle_band(levels_by_L[Lbig], E0_by_L[Lbig], Lbig)
    fit = fit_dispersion_fourier(kvals, epsvals; rmax = config.rmax_dispersion)
    println("\nε(k) Fourier coefficients (c0..c$(config.rmax_dispersion)): ",
            round.(fit.c, digits = 4))

    exact_disp = benchmark == :free_fermion ? (k -> exact_free_fermion_dispersion(params, k)) :
                 benchmark == :xxz          ? (k -> xxz_magnon_dispersion(params, k)) :
                 benchmark == :hubbard      ? (k -> hubbard_dispersion(params, k)) : nothing
    savefig(plot_dispersion(kvals, epsvals, fit; exact = exact_disp,
                            title = "[$name] one-particle dispersion (L=$Lbig)"),
            joinpath(pdir, "dispersion.png"))

    if benchmark == :free_fermion
        err = maximum(abs(epsilon(fit, k) - exact_free_fermion_dispersion(params, k))
                      for k in kvals)
        println(">>> FREE-FERMION check: max |ε_fit − ε_exact| over data = ",
                round(err, sigdigits = 3))
    end

    m1 = minimum(epsilon(fit, k) for k in range(0, π; length = 400))   # one-particle gap

    # --- (3,4) two-particle matching + phase-shift extraction -------------
    all_matches = TwoParticleMatch[]
    all_points  = PhaseShiftPoint[]
    for L in config.Ls
        for Kint in config.Ksectors
            matches = match_two_particle_levels(levels_by_L[L], E0_by_L[L], fit, L, Kint;
                                                tolerance = 0.2,
                                                mmax = min(config.max_two_particle_levels, L ÷ 2))
            append!(all_matches, matches)

            # two-particle comparison plot
            sec = levels_in_sector(levels_by_L[L], Kint)
            skip = (Kint == 0) ? 2 : 1
            numeric_ΔE = [lv.energy - E0_by_L[L] for lv in sec if lv.index > skip]
            frees = free_two_particle_energies(fit, L, Kint; mrange = 1:(L ÷ 2))
            savefig(plot_two_particle_levels(numeric_ΔE, frees, matches, Kint, L),
                    joinpath(pdir, "twoparticle_L$(L)_K$(Kint).png"))

            append!(all_points,
                    phase_points_for_sector(levels_by_L[L], E0_by_L[L], fit, L, Kint;
                                            tolerance = 0.2,
                                            mmax = min(config.max_two_particle_levels, L ÷ 2)))
        end
    end
    points_all = unwrap_phase_shifts(all_points)                      # everything (→ CSV)
    # cleanest Lüscher points: the lowest two-particle levels per (L,K)
    points = unwrap_phase_shifts(keep_lowest_per_sector(all_points, keep_lowest))

    # --- (5) benchmark overlay for the phase shift ------------------------
    bench = nothing       # δ(q) reference curve
    sbench = nothing      # S(q) reference curve (branch-free)
    golden_measured = NaN
    if benchmark == :e8
        qs = range(0.01, π; length = 200)
        _, δs = e8_phase_curve(fit, qs, m1)
        bench = (collect(qs), δs, "exact E8 δ₁₁(θ(q))")
        Sre = [real(e8_S11(rapidity_from_q(fit, q, m1))) for q in qs]
        Sim = [imag(e8_S11(rapidity_from_q(fit, q, m1))) for q in qs]
        sbench = (collect(qs), Sre, Sim, "exact E8 S₁₁")

        # golden-ratio mass ratio m2/m1, shown for every L to expose convergence
        println("\n>>> E8 golden-ratio check  m₂/m₁  (exact = 1.6180):")
        for L in config.Ls
            sec0 = levels_in_sector(levels_by_L[L], 0)
            length(sec0) ≥ 3 || continue
            e1 = sec0[2].energy - E0_by_L[L]      # m1 (A1 at rest)
            e2 = sec0[3].energy - E0_by_L[L]      # m2 (A2 at rest)
            ratio = e2 / e1
            below = e2 < 2e1 ? "" : "  (⚠ not below 2m₁ threshold — A2 unresolved)"
            println("    L=$L :  m1=", round(e1, sigdigits = 4),
                    "  m2=", round(e2, sigdigits = 4),
                    "  m2/m1=", round(ratio, digits = 4),
                    "  (err ", round(100 * abs(ratio - 1.618034) / 1.618034, sigdigits = 2), "%)", below)
            L == Lbig && (golden_measured = ratio)
        end
    elseif benchmark == :xxz
        # EXACT Bethe-ansatz two-magnon phase shift / S-matrix (continuous curve)
        qs = range(0.001, π - 0.001; length = 300)
        bench  = (collect(qs), [xxz_two_magnon_delta(params, q) for q in qs], "exact XXZ δ(q)")
        Sre = [real(xxz_two_magnon_S(params, q)) for q in qs]
        Sim = [imag(xxz_two_magnon_S(params, q)) for q in qs]
        sbench = (collect(qs), Sre, Sim, "exact XXZ S(q)")
        if !isempty(points)
            errs = [abs((p.S_real + im * p.S_imag) - xxz_two_magnon_S(params, p.q)) for p in points]
            println("\n>>> XXZ exact-S check (lowest $keep_lowest levels/L): " *
                    "max |S_extr − S_exact| = ", round(maximum(errs), sigdigits = 3),
                    " ,  mean = ", round(sum(errs)/length(errs), sigdigits = 3))
        end
    elseif benchmark == :hubbard
        # EXACT two-electron (singlet) S-matrix / phase shift
        qs = range(0.001, π - 0.001; length = 300)
        bench  = (collect(qs), [hubbard_singlet_delta(params, q) for q in qs], "exact Hubbard δ(q)")
        Sre = [real(hubbard_singlet_S(params, q)) for q in qs]
        Sim = [imag(hubbard_singlet_S(params, q)) for q in qs]
        sbench = (collect(qs), Sre, Sim, "exact Hubbard S(q)")
        if !isempty(points)
            errs = [abs((p.S_real + im * p.S_imag) - hubbard_singlet_S(params, p.q)) for p in points]
            println("\n>>> HUBBARD exact-S check (lowest $keep_lowest levels/L): " *
                    "max |S_extr − S_exact| = ", round(maximum(errs), sigdigits = 3),
                    " ,  mean = ", round(sum(errs)/length(errs), sigdigits = 3))
        end
    elseif benchmark == :free_fermion
        # exact result: S = −1 (free fermions) ⇒ δ = π/2 (mod π), independent of q
        if !isempty(points)
            qmin, qmax = extrema(p.q for p in points)
            bench = ([qmin, qmax], [π/2, π/2], "exact: S=−1 ⇒ δ=π/2 (mod π)")
            sbench = ([qmin, qmax], [-1.0, -1.0], [0.0, 0.0], "exact S=−1")
        end
    end

    # for the integrable benchmarks, fold the (mod-π) extracted δ to the branch
    # nearest the exact curve so the δ(q) plot is directly comparable
    fold_ref = benchmark == :xxz          ? (q -> xxz_two_magnon_delta(params, q)) :
               benchmark == :hubbard      ? (q -> hubbard_singlet_delta(params, q)) :
               benchmark == :free_fermion ? (q -> π/2) : nothing
    savefig(plot_phase_shift(points; benchmark = bench, fold_ref = fold_ref,
                             title = "[$name] phase shift δ_K(q)  (lowest $keep_lowest levels/L)"),
            joinpath(pdir, "phaseshift.png"))
    savefig(plot_two_delta(points; title = "[$name] 2 δ_K(q)"),
            joinpath(pdir, "twodelta.png"))
    if !isempty(points)
        savefig(plot_smatrix(points; benchmark = sbench, title = "[$name] S_K(q)"),
                joinpath(pdir, "smatrix.png"))
    end

    # --- diagnostic: ALL extracted levels (no keep_lowest filter) -------------
    # Shows the multi-particle / inelastic contamination that the clean plots
    # avoid by keeping only the lowest levels per L.
    if !isempty(points_all)
        savefig(plot_phase_shift(points_all; benchmark = bench, fold_ref = fold_ref,
                                 title = "[$name] δ_K(q) — ALL levels (contaminated)"),
                joinpath(pdir, "phaseshift_all.png"))
        savefig(plot_smatrix(points_all; benchmark = sbench,
                             title = "[$name] S_K(q) — ALL levels (contaminated)"),
                joinpath(pdir, "smatrix_all.png"))
    end

    # --- multi-K view: phase shift at SEVERAL total momenta K = 2π Kint/L -----
    # The pipeline is K-general; here we extract a few nonzero-K sectors too.
    # δ_K(q) is L-dependent for K≠0 (K=2π Kint/L), so this is done at ONE L —
    # the largest L NOT divisible by 4 (which would hit the q=π/2 artifact).
    if !isempty(multiK)
        cand = filter(L -> L % 4 != 0, config.Ls)
        Lmk  = isempty(cand) ? maximum(config.Ls) : maximum(cand)
        mkraw = PhaseShiftPoint[]
        for Kint in multiK
            append!(mkraw, phase_points_for_sector(levels_by_L[Lmk], E0_by_L[Lmk], fit, Lmk, Kint;
                              tolerance = 0.3, mmax = Lmk ÷ 2, verbose = false))
        end
        mkpts = unwrap_phase_shifts(mkraw)        # ALL elastic scattering levels per K
        de = delta_exact_K === nothing ? nothing : ((Kint, q) -> delta_exact_K(params, 2π * Kint / Lmk, q))
        fr = delta_exact_K === nothing ? π/2 : nothing       # fold near π/2 when no exact curve
        savefig(plot_phase_shift_multiK(mkpts; delta_exact = de, fold_ref = fr,
                  title = "[$name] δ_K(q) at L=$Lmk — total momenta Kint ∈ {$(join(multiK, ", "))}"),
                joinpath(pdir, "phaseshift_multiK.png"))
        println(">>> multi-K plot saved (L=$Lmk, Kint=$(multiK))")
    end

    # --- Wigner time delay  τ(q) = 2 dδ/dE :  exact derivative vs data --------
    delta_exact = benchmark == :xxz          ? (q -> xxz_two_magnon_delta(params, q)) :
                  benchmark == :hubbard      ? (q -> hubbard_singlet_delta(params, q)) :
                  benchmark == :free_fermion ? (q -> π/2) :
                  benchmark == :e8           ? (q -> e8_delta11(rapidity_from_q(fit, q, m1))) : nothing
    eps_exact   = benchmark == :xxz          ? (q -> xxz_magnon_dispersion(params, q)) :
                  benchmark == :hubbard      ? (q -> hubbard_dispersion(params, q)) :
                  benchmark == :free_fermion ? (q -> exact_free_fermion_dispersion(params, q)) :
                                               (q -> epsilon(fit, q))   # fitted ε otherwise
    if delta_exact !== nothing && !isempty(points)
        qs_td = range(0.1, π - 0.05; length = 200)
        exact_td = time_delay_exact(delta_exact, eps_exact, qs_td)
        data_td  = time_delay_points(points, fit)                  # central differences
        savefig(plot_time_delay(data_td, exact_td;
                                title = "[$name] Wigner time delay τ(q)=2 dδ/dE"),
                joinpath(pdir, "timedelay.png"))
    end

    # --- (6) global spectral fit  δ(q;θ)=a0+a1 q+a2 q²+a3 q³ ---------------
    θ = Float64[]; cost = NaN
    if length(all_matches) ≥ 2
        θ, cost = fit_phase_shift_global(all_matches, fit; order = 3)
        println("\nGlobal phase-shift fit  δ(q)=a0+a1 q+a2 q²+a3 q³")
        println("    θ = ", round.(θ, sigdigits = 4), "   cost D = ", round(cost, sigdigits = 4))
        savefig(plot_global_fit(all_matches, fit, θ;
                                title = "[$name] global fit: numerical vs predicted"),
                joinpath(pdir, "globalfit.png"))

        # 1D scan of the constant phase a0, and 2D scan of (a0,a1)
        a0s = range(-π, π; length = 121)
        v, negD = scan_delta_1d(all_matches, fit, 1, a0s; base = θ)
        savefig(plot_scan_1d(v, negD; label = "a0", title = "[$name] scan −D(a0)"),
                joinpath(pdir, "scan_a0.png"))
        a1s = range(-2.0, 2.0; length = 81)
        vi, vj, G = scan_delta_2d(all_matches, fit, (1, 2), a0s, a1s; base = θ)
        savefig(plot_scan_2d(vi, vj, G; xlabel = "a0", ylabel = "a1",
                             title = "[$name] scan −D(a0,a1)"),
                joinpath(pdir, "scan_a0a1.png"))
    end

    # --- (7) save results -------------------------------------------------
    write_phase_csv(joinpath(rdir, "phase_points.csv"), points_all)
    summary = Dict{String, Any}(
        "scenario"        => name,
        "model"           => string(typeof(params)),
        "params"          => string(params),
        "Ls"              => Float64.(config.Ls),
        "Ksectors"        => Float64.(config.Ksectors),
        "rmax_dispersion" => Float64(config.rmax_dispersion),
        "eps_fit_coeffs"  => fit.c,
        "m1_gap"          => m1,
        "n_phase_points"  => Float64(length(points_all)),
        "global_fit_theta" => isempty(θ) ? Float64[] : θ,
        "global_fit_cost" => cost,
        "golden_ratio_measured" => golden_measured,
        "golden_ratio_exact"    => 1.6180339887,
    )
    write_json(joinpath(rdir, "summary.json"), summary)
    println("\nSaved: plots → examples/luscher/plots/$(name)/ ; " *
            "results → examples/luscher/results/$(name)/")
    return (; fit, points, matches = all_matches, levels_by_L, E0_by_L, m1, theta = θ)
end

# ============================================================================
#  Run the scenarios
# ============================================================================
results = Dict{String, Any}()

# (A) XXZ easy-axis ferromagnet — HEADLINE analytic benchmark.
#     Bethe-ansatz integrable: the extracted δ(q)/S(q) are compared against the
#     EXACT closed-form two-magnon S-matrix S(q) = −(1−Δe^{−iq})/(1−Δe^{iq}).
results["xxz"] = analyze("xxz",
    XXZParams(J = 1.0, Δ = 1.5, h = 2.0),
    LuscherConfig(Ls = [10, 12, 14, 16], Ksectors = [0],
                  nlevels = 16, rmax_dispersion = 5, max_two_particle_levels = 8);
    benchmark = :xxz, keep_lowest = 2,   # 2 lowest scattering levels/L: clean q∈[0.4,1.3]
                                         # (avoids the q=π/2 commensurate-lattice artifact at L%4==0)
    multiK = [0, 1, 2, 3],               # nonzero total momentum: δ_K(q) vs EXACT K-dependent
    delta_exact_K = xxz_two_magnon_delta)  # Bethe S-matrix S_K(q)=−(cos(K/2)−Δe^{−iq})/(cos(K/2)−Δe^{iq})

# (A2) Hubbard two-electron benchmark — exact Lieb–Wu singlet S-matrix.
#      Empty vacuum; isolate the 1↑1↓ singlet via Sᶻ/N penalties + K=0 reflection
#      parity.  d=4 ⇒ ED caps at small L (we use L=6,7,8,9).
results["hubbard"] = analyze("hubbard",
    HubbardParams(t = 1.0, U = 4.0, μ = 3.0, Λ = 10.0),
    LuscherConfig(Ls = [6, 7, 8, 9], Ksectors = [0],
                  nlevels = 12, rmax_dispersion = 4, max_two_particle_levels = 6);
    benchmark = :hubbard, keep_lowest = 2)

# (B) Ising free-fermion validation (S = −1, exact dispersion).
#     Deep in the disordered phase (h_x=4 ≫ h_c=1): gap 2(h_x−1)=6 exceeds the
#     bandwidth ≈4, so the one-particle band is well separated from the 4-particle
#     continuum and the low two-particle levels are clean δ=π/2 (near the critical
#     point h_x≈1.5 the bands overlap and multi-particle states contaminate).
results["freefermion"] = analyze("freefermion",
    IsingParams(J = 1.0, hx = 4.0, hz = 0.0),
    LuscherConfig(Ls = [8, 10, 12, 14], Ksectors = [0],
                  nlevels = 14, rmax_dispersion = 4, max_two_particle_levels = 6);
    benchmark = :free_fermion,
    multiK = [0, 1, 2],                          # S=−1 ⇒ δ=π/2 at every total momentum
    delta_exact_K = (p, K, q) -> π/2)

# (B2) Ising DISORDERED + longitudinal field — interacting, NON-integrable, but
#      with a cleanly SEPARATED one-particle band (unlike the critical E8 point).
#      Same h_x=4 as the free-fermion point (gap≈6.8 > bandwidth≈3.6, band fully
#      below the 2-particle threshold), now with h_z=1 switched on: the magnons
#      genuinely interact, so the extracted δ_K(q) deviates from the free π/2 in a
#      smooth, clean way.  No closed-form S-matrix (non-integrable), so we fold the
#      multi-K curves near π/2 to expose the K-dependent interaction shift.
results["disordered"] = analyze("disordered",
    IsingParams(J = 1.0, hx = 4.0, hz = 1.0),
    LuscherConfig(Ls = [8, 10, 12, 14], Ksectors = [0],
                  nlevels = 14, rmax_dispersion = 4, max_two_particle_levels = 6);
    benchmark = :none, keep_lowest = 4,
    multiK = [0, 1, 2, 3])               # K-resolved δ even without an exact S-matrix

# (C) Ising E8 point — literature benchmark: golden-ratio mass ratio m₂/m₁
results["e8"] = analyze("e8",
    IsingParams(J = 1.0, hx = 1.0, hz = 0.04),
    LuscherConfig(Ls = [8, 10, 12, 14], Ksectors = [0],
                  nlevels = 16, rmax_dispersion = 4, max_two_particle_levels = 6);
    benchmark = :e8, keep_lowest = 1)

# (C2) Ising E8 with a LARGER longitudinal field — isolates the one-particle band.
#      The E8 *field theory* is integrable for any magnetic coupling, but on the
#      LATTICE E8 is the small-h_z scaling limit.  Raising h_z to 0.5 pushes the band
#      top below the two-particle threshold (gap m₁≈3.7, band top≈5.4 < 2m₁≈7.4), so
#      the PHASE extraction is clean — at the cost of leaving the precise scaling
#      limit: the golden ratio drifts 0.09%→~1.4% and the phase sits ~0.1–1 rad off
#      E8 (a crossover, "approximately E8 with an isolated band", not exact E8).
results["e8_isolated"] = analyze("e8_isolated",
    IsingParams(J = 1.0, hx = 1.0, hz = 0.5),
    LuscherConfig(Ls = [10, 12, 14, 16], Ksectors = [0],
                  nlevels = 16, rmax_dispersion = 4, max_two_particle_levels = 6);
    benchmark = :e8, keep_lowest = 2)

# (D) Ising assignment defaults — non-integrable confinement (physics, no exact δ)
results["confinement"] = analyze("confinement",
    IsingParams(J = 1.0, hx = 0.6, hz = 0.05),
    LuscherConfig(Ls = [10, 12], Ksectors = [0],
                  nlevels = 12, rmax_dispersion = 4, max_two_particle_levels = 5);
    benchmark = :none)

println("\nAll scenarios complete.")
