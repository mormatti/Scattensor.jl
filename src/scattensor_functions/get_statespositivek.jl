"""
    get_statespositivek(states::Vector{<:BlochState}) -> Vector{<:BlochState}

Filter Bloch states by momentum, keeping those with nonnegative momentum `k ≥ 0`.

# Arguments
- `states`: Vector of Bloch states.

# Returns
- A vector containing only states with `momentum(state) ≥ 0`.

# Notes
- This uses `momentum(state)`, i.e. `π * state.koverpi`. No Brillouin-zone wrapping is applied.
"""
function get_statespositivek(states::Vector{<:BlochState})
    return filter(state -> momentum(state) >= 0, states)
end
