# TODO: find a way to write these function in a more efficient 
# (complexity) and compact way

function get_groundstate(states::Vector{BlochStateType}) where {BlochStateType <: BlochState}

    # We identify the position of the state which have the lowest energy
    i0 = 1
    for i in eachindex(states)
        if energy(states[i]) < energy(states[i0])
            i0 = i
        end
    end
    return states[i0]
end
export get_groundstate