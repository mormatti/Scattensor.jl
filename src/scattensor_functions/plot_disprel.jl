function plot_disprel(disprelvec::Vector{T}; kwargs...) where {T <: BlochState}
    E = [energy(state) for state in disprelvec]
    k = [momentum(state) for state in disprelvec]
    return Plots.scatter(k, E; kwargs...)
end

export plot_disprel