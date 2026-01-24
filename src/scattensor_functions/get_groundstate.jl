"""
    get_groundstate(states::Vector{<:BlochState}) -> BlochState

Return the element of `states` with the lowest energy.

# Arguments
- `states`: Vector of `BlochState` objects (or subtypes).

# Returns
- The state `s ∈ states` that minimizes `energy(s)`.

# Notes
- This is an `O(N)` scan over `states`.
- If `states` is empty, this function throws a bounds error (consistent with the current
  implementation). Consider guarding with `isempty(states)` in user code if needed.
"""
function get_groundstate(states::Vector{<:BlochState})
    igs = 1
    for i in eachindex(states)
        if energy(states[i]) < energy(states[igs])
            igs = i
        end
    end
    return states[igs]
end

export get_groundstate