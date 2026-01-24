"""
    tdvp_time_evolution!(data, H, ψ, dt, Δt, H0; cutoff=default_cutoff, maxdim=default_maxdim, groundstate=:none, kwargs...) -> nothing

Time-evolve an MPS with TDVP and record diagnostics into `data`.

At each step, this routine applies `ITensorMPS.tdvp` to evolve the state by `-im*dt`, normalizes the
result, stores intermediate states and diagnostics, and writes a small set of `.png` plots to disk
(times, estimated remaining times, memory size, link dimensions, and local energies).

# Arguments
- `data::Dict`: Output dictionary mutated in-place to store results and metadata.
- `H::MPO`: Time-evolution Hamiltonian MPO (must have same length as `ψ`).
- `ψ::MPS`: Initial state (modified by evolution, but the evolved states are stored in `data[:states]`).
- `dt::Real`: Time step.
- `Δt::Real`: Total evolution time (number of steps is `round(Int, Δt/dt)`).
- `H0::MPO`: Local Hamiltonian used to compute local energy densities via `local_expvals`.

# Keyword Arguments
- `cutoff=default_cutoff`: TDVP truncation cutoff.
- `maxdim=default_maxdim`: Maximum bond dimension during TDVP.
- `groundstate=:none`: If not `:none`, interpreted as an MPS used to compute a reference local energy density
  (subtracted from the energies measured during evolution).
- `kwargs...`: Forwarded to `tdvp`.

# Returns
- `nothing`. Results are stored in `data`.

# Side effects
- Writes several PNG files in the current working directory.

# Notes
- This function is marked TODO in the source (“tidy and document”) and should be considered experimental.
"""

function tdvp_time_evolution!(data::Dict, H::MPO, ψ::MPS, dt::Real, Δt::Real, H0::MPO; cutoff = default_cutoff, maxdim = default_maxdim, groundstate = :none, kwargs...)
    
    # Change the backend to Plots.jl
    gr()

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

    if groundstate != :none
        println("Computing ground state for energy reference...")
        gsenergies = local_expvals(groundstate, H0, hermitian = true)
    else
        gsenergies = zeros(L - size(H)[1] + 1)
    end

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
            ψ = tdvp(H, -im * dt, ψ; cutoff = cutoff, maxdim = maxdim, kwargs...)
            normalize!(ψ)
            push!(data[:states], ψ)
            push!(data[:linkdimens], linkdims(ψ))
            push!(data[:energies], local_expvals(ψ, H0, hermitian = true) .- gsenergies)
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
