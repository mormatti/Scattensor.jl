function partial_trace(mpo::MPO, j1::Integer, j2::Integer)
    @assert j1 > 0 && j1 < length(mpo) + 1 "j1 out of range"
    @assert j2 > 0 && j2 < length(mpo) + 1 "j2 out of range"
    @assert j1 < j2 "j1 must be less than j2"
    if j1 == 1 && j2 == length(mpo)
        return mpo
    end
    if j1 == 1
        return partial_trace_right(mpo, j2)
    end
    if j2 == length(mpo)
        return partial_trace_left(mpo, j1)
    else
        return partial_trace_left(partial_trace_right(mpo, j2), j1)
    end
end

function _partial_trace_left(mpo::MPO, j1::Integer)
    @assert j1 > 0 && j1 < length(mpo) + 1 "j1 out of range"
    if j1 == 1
        return mpo
    end
    sites = siteinds(mpo)
    A = mpo[1] * delta(sites[1][1], sites[1][2])
    for j in 2:j1-1
        A = A * mpo[j] * delta(sites[j][1], sites[j][2])
    end
    # We consider the MPO M from index j1 to the end
    tensorlist = deepcopy(mpo[j1:end])
    tensorlist[1] = A * mpo[j1]
    return MPO(tensorlist)
end

function _partial_trace_right(mpo::MPO, j1::Integer)
    @assert j1 > 0 && j1 < length(mpo) + 1 "j1 out of range"
    if j1 == length(mpo)
        return mpo
    end
    sites = siteinds(mpo)
    A = mpo[end] * delta(sites[end][1], sites[end][2])
    for j in length(mpo)-1:-1:j1+1
        A = mpo[j] * delta(sites[j][1], sites[j][2]) * A
    end
    # We consider the MPO M from the beginning to index j1
    tensorlist = deepcopy(mpo[1:j1])
    tensorlist[end] = mpo[j1] * A
    return MPO(tensorlist)
end

export partial_trace