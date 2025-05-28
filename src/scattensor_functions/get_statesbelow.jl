function get_statesbelow(states::Vector{BlochStateType}, enrg::RealType) where {BlochStateType <: BlochState, RealType <: Real}
    return filter(state -> energy(state) <= enrg, states)
end

export get_statesbelow