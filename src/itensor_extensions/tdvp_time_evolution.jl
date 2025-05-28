function tdvp_time_evolution!(
    data::Dict, 
    H::MPO,
    ψ::MPS, 
    dt::dtType, 
    Δt::ΔtType, 
    H0::MPO; 
    kwargs...) where {
        dtType <: Real, 
        ΔtType <: Real
        }
    
    # Getting the local dimension, ensuring uniformity
    sites = siteinds(ψ)
    site_dims = [dim(sites[i]) for i in eachindex(sites)]
    all_same = length(unique(site_dims)) == 1
    if all_same
        data[:d] = site_dims[1]
    else
        # TODO: we should allow different dimensions actually...
        error("The site dimensions are not all the same.")
    end

    data[:states] = Vector{MPS}()
    @assert length(H) == length(ψ)
    data[:dt] = dt
    data[:Δt] = Δt
    data[:kwargs] = kwargs
    L = length(H)
    data[:L] = L
    N = Integer(round(Δt / dt))
    data[:N] = Integer(round(Δt / dt))
    simulationtime = 0
    data[:simulationtime] = 0
    times = [n * dt for n in 1:N]
    data[:times] = times
    positions = [x for x in 2:(L-1)]
    data[:positions] = positions

    println("Starting the simulation...")
    println("Total number of steps to execute: $N")

    # Matrices
    data[:linkdimens] = []
    data[:energies] = []
    # Arrays
    data[:timeselapsed] = []
    data[:estimatedtimes] = []
    data[:totsizes] = []

    for n in 1:N

        timelapsed = @elapsed begin
            ψ = tdvp(H, -im * dt, ψ; kwargs...)
            normalize!(ψ)

            push!(data[:states], ψ)
            push!(data[:linkdimens], linkdims(ψ))
            push!(data[:energies], local_expvals(ψ, H0, data[:d]))
        end
        push!(data[:timeselapsed], timelapsed)
        push!(data[:estimatedtimes], timelapsed * (N - n))
        push!(data[:totsizes], Base.summarysize(data[:states]))

        pl1 = Plots.plot(data[:timeselapsed], dpi=300)
        Plots.savefig(pl1, "times_elapsed.png")
        pl2 = Plots.plot(data[:estimatedtimes], dpi=300)
        Plots.savefig(pl2, "estimated_times.png")
        pl3 = Plots.plot(data[:totsizes], dpi=300)
        Plots.savefig(pl3, "total_sizes.png")
        pl4 = Plots.heatmap(hcat(data[:linkdimens]...), dpi=300)
        Plots.savefig(pl4, "link_dimensions.png")
        pl5 = Plots.heatmap(hcat(data[:energies]...), dpi=300)
        Plots.savefig(pl5, "local_energies.png")
    end

    data[:simulationtime] = simulationtime
end

export tdvp_time_evolution!