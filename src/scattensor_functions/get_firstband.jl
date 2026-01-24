"""
    get_firstband(states::Vector{<:BlochState}) -> Vector{<:BlochState}

Construct a "first band" by selecting, for each distinct momentum `k`, the lowest-energy state.

This is a convenience routine for post-processing the output of `dispersion_relation`, which often
returns multiple eigenstates at each momentum.

# Arguments
- `states`: Vector of Bloch states (typically containing multiple eigenstates per momentum).

# Returns
- A vector containing exactly one state per momentum: the minimum-energy one.
  The output is ordered by increasing momentum.

# Notes
- Momentum grouping is performed by equality comparisons of `momentum(state)`. For robust grouping,
  prefer storing `koverpi` as a `Rational` in your `BlochState`s (floating-point momentum grids may
  suffer from equality/roundoff issues).
"""
function get_firstband(states::Vector{BlochStateType}) where BlochStateType <: BlochState
    # If the list of states is empty, we return an empty list    
    if isempty(states)
        return Vector{BlochStateType}()
    else
        # We get the list of all the differents ks inside states
        ks = unique(momentum.(states))
        # We sort the list of ks
        sort!(ks)
        # For each k, we get the state with the lowest energy
        # The following 2 lines are important to initialize band as a Vector of the correct type
        band = [get_groundstate(states)]
        pop!(band)
        for k in ks
            states_k = filter(state -> momentum(state) == k, states)
            ground_k = get_groundstate(states_k)
            push!(band, ground_k)
        end
        return band
    end
end

export get_firstband