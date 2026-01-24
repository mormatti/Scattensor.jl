"""
    get_statesabove(states::Vector{<:BlochState}, enrg::Real) -> Vector{<:BlochState}

Filter Bloch states by energy threshold, keeping those with `energy(state) ≥ enrg`.

# Arguments
- `states`: Vector of Bloch states.
- `enrg`: Energy cutoff.

# Returns
- A vector of all states in `states` whose energy is greater than or equal to `enrg`.
"""
function get_statesabove(states::Vector{<:BlochState}, enrg::Real)
    return filter(state -> energy(state) >= enrg, states)
end

export get_statesabove