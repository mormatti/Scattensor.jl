"""
    pop_groundstate!(states::Vector{<:BlochState}) -> BlochState

Remove and return the lowest-energy state from `states`.

# Arguments
- `states`: Vector of `BlochState` objects (or subtypes). The vector is modified in-place.

# Returns
- The state with minimal `energy` among those originally in `states`.

# Notes
- This is an `O(N)` scan over `states` to locate the minimum-energy entry.
- If `states` is empty, this function throws a bounds error (consistent with the current implementation).
"""
function pop_groundstate!(states::Vector{<:BlochState})
    # We identify the position of the state which have the lowest energy
    i0 = 1
    for i in eachindex(states)
        if energy(states[i]) < energy(states[i0])
            i0 = i
        end
    end
    return popat!(states, i0)
end

export pop_groundstate!