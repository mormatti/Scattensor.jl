function Plots.plot(disprelvec::Vector{<:BlochState}; kwargs...)
    E = [energy(state) for state in disprelvec]
    k = [momentum(state) for state in disprelvec]
    return Plots.scatter(k, E; kwargs...)
end

export plot_disprel