function get_statesabove(states::Vector{<:BlochState}, enrg::Real)
    return filter(state -> energy(state) >= enrg, states)
end

export get_statesabove