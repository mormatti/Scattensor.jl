# TODO: write a documentation for this function.
function pop_firstband!(states::Vector{BlochStateType}) where BlochStateType <: BlochState
    # If the list of states is empty, we return an empty list
    if isempty(states)
        return Vector{BlochStateType}()
    else
        # We get the list of all the differents ks inside states
        ks = unique(momentum.(states))
        # We sort the list of ks
        sort!(ks)
        # For each k, we get the state with the lowest energy, and we remove it from states
        # The following 2 lines are important to initialize band as a Vector of the correct type
        # Actually, you sould fix that.
        band = [get_groundstate(states)]
        pop!(band)
        for k in ks
            states_k = filter(state -> momentum(state) == k, states)
            ground_k = pop_groundstate!(states_k)
            push!(band, ground_k)
            # We remove the state from states
            deleteat!(states, findfirst(==(ground_k), states))
        end
        return band
    end
end

export pop_firstband!