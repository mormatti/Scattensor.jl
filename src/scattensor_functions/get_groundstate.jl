# TODO: find a way to write these function in a more efficient (complexity) and compact way
# TODO: write a documentation for this function.
function get_groundstate(states::Vector{<:BlochState})
    igs = 1
    for i in eachindex(states)
        if energy(states[i]) < energy(states[igs])
            igs = i
        end
    end
    return states[igs]
end

export get_groundstate