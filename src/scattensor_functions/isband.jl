# We write a function that given a list of states it corresponds to a band.
# Namely, it checks if for each k corresponds to only one state.
function isband(states::Vector{BlochStateType}) where BlochStateType <: BlochState
    ks = momentum.(states)
    unique_ks = unique(ks)
    if length(ks) != length(unique_ks)
        return false
    else
        return true
    end
end

export isband