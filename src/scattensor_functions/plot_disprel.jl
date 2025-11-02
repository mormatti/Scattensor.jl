using Plots
using LaTeXStrings

function plot_disprel(disprelvec::Vector{<:BlochState}, from::Real, to::Real; kwargs...)

    Elist = []
    klist = []

    function catch_ks(koverpi, from, to)
        koverpicopy = koverpi
        koverpis = []
        while koverpicopy * π >= from
            koverpicopy -= 2
        end
        while koverpicopy * π <= to
            if from <= koverpicopy * π <= to
                push!(koverpis, koverpicopy)
            end
            koverpicopy += 2
        end
        return koverpis
    end

    for state in disprelvec
        for kvp in catch_ks(state.koverpi, from, to)
            push!(klist, kvp * π)
            push!(Elist, energy(state))
        end
        for kvp in catch_ks(-state.koverpi, from, to)
            push!(klist, kvp * π)
            push!(Elist, energy(state))
        end
        # TODO for parity invariance please complete here
    end

    Emin, Emax = extrema(Elist)

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
    klims = padded_limits(from, to)
    Elims = padded_limits(Emin, Emax)

    # Canonical tick candidates and labels
    tick_pos_all  = [-2π, -(3/2) * π, -π, -π/2, 0, π/2, π, (3/2) * π, 2π]
    tick_lab_all  = [L"2\pi", L"-3\pi/2", L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"]

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

    return Plots.scatter(klist, Elist; defaults..., kwargs...)
end

export plot_disprel