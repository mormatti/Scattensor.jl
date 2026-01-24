using ITensors

"""
    CustomObserver(tolerance)

Simple DMRG observer that stops when energy converges.

This `AbstractObserver` implementation can be passed to `ITensorMPS.dmrg` and will request early
termination when the absolute change in energy between sweeps drops below `tolerance`.

# Arguments
- `tolerance`: Convergence threshold on `abs(energy - lastenergy)`.

# Notes
- This observer also prints a small progress message on the first sweep.
"""
mutable struct CustomObserver <: AbstractObserver
    lastenergy::Float64
    tolerance::Float64
end

CustomObserver(tolerance) = CustomObserver(Inf, tolerance)

# The following function is called at the end of each sweep
function ITensorMPS.checkdone!(observer::CustomObserver; energy, sweep, kwargs...)

    if sweep == 1
        println("Computing DMRG.")
    end

    print("\r", "\x1b[2A", "\x1b[1M", "\n")

    stop = (sweep > 1) && (abs(energy - observer.lastenergy) < observer.tolerance)
  
    if stop
        # Here print something at the end
        # println("Final energy = $(abs(energy - observer.lastenergy)) < $(observer.tolerance): stopping after sweep $sweep")
    end

    observer.lastenergy = energy

    return stop
end