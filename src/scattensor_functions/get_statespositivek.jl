# TODO: write a documentation for this function.
function get_statespositivek(states::Vector{<:BlochState})
    return filter(state -> momentum(state) >= 0, states)
end

export get_statespositivek