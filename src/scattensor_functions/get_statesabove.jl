# TODO: write a documentation for this function.
function get_statesabove(states::Vector{<:BlochState}, enrg::Real)
    return filter(state -> energy(state) >= enrg, states)
end

export get_statesabove