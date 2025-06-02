function get_firstband(states::Vector{BlochStateType}) where BlochStateType <: BlochState    
    if isempty(states)
        return Vector{BlochStateType}()
    else
        band = [get_groundstate(states)]
        pop!(band)
        while !isempty(states)
            ground = get_groundstate(states)
            k = momentum(ground)
            # We remove all the states which have momentum k
            states = filter(state -> momentum(state) != k, states)
            push!(band, ground)
        end
        return band
    end
end

export get_firstband