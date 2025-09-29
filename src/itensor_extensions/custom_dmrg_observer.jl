using ITensors

mutable struct CustomObserver <: AbstractObserver
    lastenergy::Float64
    tolerance::Float64
end

CustomObserver(tolerance) = CustomObserver(Inf, tolerance)

# The following function is called at the end of each sweep
function ITensorMPS.checkdone!(observer::CustomObserver; energy, sweep, outputlevel = 1, kwargs...)

    if sweep == 1
        println("Computing DMRG.")
    end

    cancel_terminal_line()

    println("Energy = $energy, sweep = $sweep")

    stop = (sweep > 1) && (abs(energy - observer.lastenergy) < observer.tolerance)
  
    if outputlevel > 0 && stop
        println("Computation done!")
        println("Final energy = $(abs(energy - observer.lastenergy)) < $(observer.tolerance): stopping after sweep $sweep")
    end

    observer.lastenergy = energy

    return stop
end