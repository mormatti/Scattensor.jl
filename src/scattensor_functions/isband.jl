"""
    isband(states::Vector{<:BlochState}) -> Bool

Return `true` if `states` contains at most one state for each momentum `k`.

This is a lightweight structural check often used to validate that a set of Bloch states already
represents a single band (one state per momentum point).

# Arguments
- `states`: Vector of Bloch states.

# Returns
- `true` if all `momentum.(states)` values are unique, otherwise `false`.
"""
function isband(states::Vector{BlochStateType}) where BlochStateType <: BlochState
    ks = momentum.(states)
    unique_ks = unique(ks)
    if length(ks) != length(unique_ks)
        return false
    else
        return true
    end
end
