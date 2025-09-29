using Plots
using LaTeXStrings

function Plots.plot(disprelvec::Vector{<:BlochState}; kwargs...)
    E = [energy(state) for state in disprelvec]
    k = [momentum(state) for state in disprelvec]

    # Data ranges
    kmin, kmax = extrema(k)
    Emin, Emax = extrema(E)

    # Handle degenerate ranges and add ~10% padding
    function padded_limits(a, b; frac = 0.10)
        if a == b
            δ = max(abs(a), 1.0) * frac
            return (a - δ, b + δ)
        else
            m = frac * (b - a)
            return (a - m, b + m)
        end
    end
    klims = padded_limits(kmin, kmax)
    Elims = padded_limits(Emin, Emax)

    # Canonical tick candidates and labels
    tick_pos_all  = [-π, -π/2, 0, π/2, π]
    tick_lab_all  = [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]

    # Keep only ticks inside the padded limits (with a tiny tolerance)
    tol = 1e-9 * (abs(klims[2] - klims[1]) + 1)
    mask = map(x -> (klims[1] - tol) <= x <= (klims[2] + tol), tick_pos_all)
    tick_pos = [x for (x, m) in zip(tick_pos_all, mask) if m]
    tick_lab = [s for (s, m) in zip(tick_lab_all, mask) if m]

    # If none of the canonical ticks are inside, let Plots choose automatically
    xticks_setting = isempty(tick_pos) ? :auto : (tick_pos, tick_lab)

    defaults = (
        xlabel = L"k",
        ylabel = L"\mathcal{E}_k",
        aspect_ratio = 2,
        legend = false,
        xlims = klims,
        ylims = Elims,
        xticks = xticks_setting,
        framestyle = :box,
    )

    return Plots.scatter(k, E; defaults..., kwargs...)
end