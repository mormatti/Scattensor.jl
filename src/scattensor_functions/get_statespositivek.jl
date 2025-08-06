# Be careful, the momentum must be a fraction possibly... (otherwise, fix it)

function get_statespositivek(states::Vector{<:BlochState}, moment::Real)
    return filter(state -> momentum(state) >= 0, states)
end

export get_statesbelow