"""
    tdvp_time_evolution!(data, H, ψ, dt, Δt, H0;
                         cutoff = default_cutoff, maxdim = default_maxdim,
                         groundstate = :none,
                         save_states = true,
                         plots_every = 5, dashboard = true,
                         checkpoint_every = 10,
                         checkpoint_path = "tdvp_checkpoint.bin",
                         kwargs...) -> nothing

Time-evolve an MPS with TDVP and record diagnostics into `data`.

At each step the state is evolved by `-im*dt` with `ITensorMPS.tdvp`,
normalized, and measured: local energy densities (via `H0` and, if given,
referenced to `groundstate`), link dimensions, timings.

# Arguments
- `data::Dict`: output dictionary, mutated in place (keys: `:times`,
  `:energies`, `:linkdimens`, `:timeselapsed`, `:estimatedtimes`,
  `:totsizes`, and `:states` if `save_states`).
- `H::MPO`: evolution Hamiltonian. `ψ::MPS`: initial state.
- `dt`, `Δt`: time step and total time (`N = round(Δt/dt)` steps).
- `H0::MPO`: local Hamiltonian block for the energy-density measurement.

# Keyword Arguments
- `groundstate`: reference MPS whose local energies are subtracted.
- `save_states = true`: keep every evolved MPS in `data[:states]`.
  **Memory grows linearly** — set `false` for long runs (the current state is
  always available in the checkpoint).
- `plots_every = 5`: write the live PNG diagnostics every so many steps
  (0 disables plotting).
- `dashboard = true`: write `dashboard.html` in the working directory — an
  auto-refreshing page collecting the live plots (open it in a browser).
- `checkpoint_every = 10`: every so many steps, atomically serialize a
  checkpoint (0 disables). On a killed run nothing is lost but the tail.
- `checkpoint_path = "tdvp_checkpoint.bin"`: checkpoint location; read it
  back with [`tdvp_checkpoint`](@ref).
- `kwargs...`: forwarded to `ITensorMPS.tdvp`.

# Checkpoint contents
A `Dict` with `:step`, `:times`, `:energies` (matrix, one column per step),
`:linkdimens`, `:psi` (the current MPS), `:dt`, `:Δt` — enough to analyze a
partial run or restart the evolution from the last saved state.
"""
function tdvp_time_evolution!(data::Dict, H::MPO, ψ::MPS, dt::Real, Δt::Real, H0::MPO;
                              cutoff = default_cutoff, maxdim = default_maxdim,
                              groundstate = :none,
                              save_states::Bool = true,
                              plots_every::Integer = 5,
                              dashboard::Bool = true,
                              checkpoint_every::Integer = 10,
                              checkpoint_path::AbstractString = "tdvp_checkpoint.bin",
                              kwargs...)
    gr()

    sites = siteinds(ψ)
    site_dims = [dim(sites[i]) for i in eachindex(sites)]
    length(unique(site_dims)) == 1 || error("The site dimensions are not all the same.")
    data[:d] = site_dims[1]

    @assert length(H) == length(ψ)
    L = length(H)
    N = Integer(round(Δt / dt))
    data[:dt] = dt; data[:Δt] = Δt; data[:kwargs] = kwargs
    data[:L] = L; data[:N] = N
    times = [n * dt for n in 1:N]
    data[:times] = times

    gsenergies = if groundstate != :none
        println("Computing ground state for energy reference...")
        local_expvals(groundstate, H0, hermitian = true)
    else
        zeros(L - size(H)[1] + 1)
    end

    save_states && (data[:states] = Vector{MPS}())
    data[:linkdimens] = []
    data[:energies] = []
    data[:timeselapsed] = []
    data[:estimatedtimes] = []
    data[:totsizes] = []

    dashboard && _write_dashboard()

    println("Starting the simulation ($N steps)...")
    for n in 1:N
        timelapsed = @elapsed begin
            ψ = tdvp(H, -im * dt, ψ; cutoff = cutoff, maxdim = maxdim, kwargs...)
            normalize!(ψ)
            save_states && push!(data[:states], ψ)
            push!(data[:linkdimens], linkdims(ψ))
            push!(data[:energies], local_expvals(ψ, H0, hermitian = true) .- gsenergies)
        end
        push!(data[:timeselapsed], timelapsed)
        push!(data[:estimatedtimes], timelapsed * (N - n))
        push!(data[:totsizes], save_states ? Base.summarysize(data[:states]) : Base.summarysize(ψ))

        if plots_every > 0 && (n % plots_every == 0 || n == N)
            _write_live_plots(data)
        end
        if checkpoint_every > 0 && (n % checkpoint_every == 0 || n == N)
            _write_checkpoint(checkpoint_path, data, ψ, n)
        end
        print("\rstep $n / $N  (", round(timelapsed, digits = 2), " s)   ")
    end
    println()
    nothing
end

"""
    tdvp_checkpoint(path = "tdvp_checkpoint.bin") -> Dict

Load a checkpoint written by [`tdvp_time_evolution!`](@ref). Contains the
partial `:times`/`:energies`/`:linkdimens`, the step index `:step`, and the
last evolved state `:psi` (restart: pass `:psi` as the new initial state and
shorten `Δt` accordingly).
"""
tdvp_checkpoint(path::AbstractString = "tdvp_checkpoint.bin") = deserialize(path)

# ---------------------------------------------------------------------------
function _write_checkpoint(path, data, ψ, n)
    ck = Dict{Symbol, Any}(
        :step => n,
        :times => data[:times][1:n],
        :energies => hcat(data[:energies]...),
        :linkdimens => hcat(data[:linkdimens]...),
        :psi => ψ,
        :dt => data[:dt], :Δt => data[:Δt], :L => data[:L],
    )
    tmp = path * ".tmp"
    serialize(tmp, ck)
    mv(tmp, path; force = true)      # atomic on the same filesystem
end

function _write_live_plots(data)
    pl1 = Plots.plot(data[:timeselapsed], dpi = 150, legend = false, title = "seconds per step")
    Plots.savefig(pl1, "times_elapsed.png")
    pl2 = Plots.plot(data[:estimatedtimes], dpi = 150, legend = false, title = "estimated seconds remaining")
    Plots.savefig(pl2, "estimated_times.png")
    pl3 = Plots.plot(data[:totsizes], dpi = 150, legend = false, title = "memory (bytes)")
    Plots.savefig(pl3, "total_sizes.png")
    pl4 = Plots.heatmap(hcat(data[:linkdimens]...), dpi = 150, title = "link dimensions")
    Plots.savefig(pl4, "link_dimensions.png")
    pl5 = Plots.heatmap(hcat(data[:energies]...), dpi = 150, title = "local energies")
    Plots.savefig(pl5, "local_energies.png")
end

function _write_dashboard()
    html = """
    <!doctype html><html><head><meta http-equiv="refresh" content="5">
    <title>TDVP live</title></head>
    <body style="font-family:sans-serif;background:#141414;color:#eee;text-align:center">
    <h2>TDVP live dashboard <small style="color:#888">(auto-refresh 5 s)</small></h2>
    <img src="local_energies.png"  style="width:49%"><img src="link_dimensions.png" style="width:49%"><br>
    <img src="times_elapsed.png"   style="width:32%"><img src="estimated_times.png" style="width:32%"><img src="total_sizes.png" style="width:32%">
    </body></html>
    """
    write("dashboard.html", html)
    println("Live dashboard: open $(joinpath(pwd(), "dashboard.html")) in a browser.")
end
