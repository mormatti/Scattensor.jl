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