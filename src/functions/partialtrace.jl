function partial_trace_left(M::MPO, j1::j1Type) where {j1Type <: Integer}
    @assert j1 > 0 && j1 < length(M) + 1 "j1 out of range"
    if j1 == 1
        return M
    end
    sts = siteinds(M)
    A = M[1] * delta(sts[1][1], sts[1][2])
    for j in 2:j1-1
        A = A * M[j] * delta(sts[j][1], sts[j][2])
    end
    # We consider the MPO M from index j1 to the end
    M′ = deepcopy(M[j1:end])
    M′[1] = A * M[j1]
    return MPO(M′)
end

function partial_trace_right(M::MPO, j1::j1Type) where {j1Type <: Integer}
    @assert j1 > 0 && j1 < length(M) + 1 "j1 out of range"
    if j1 == length(M)
        return M
    end
    sts = siteinds(M)
    A = M[end] * delta(sts[end][1], sts[end][2])
    for j in length(M)-1:-1:j1+1
        A = M[j] * delta(sts[j][1], sts[j][2]) * A
    end
    # We consider the MPO M from the beginning to index j1
    M′ = deepcopy(M[1:j1])
    M′[end] = M[j1] * A
    return MPO(M′)
end

function partial_trace(M::MPO, j1::j1Type, j2::j2Type) where {j1Type <: Integer, j2Type <: Integer}
    @assert j1 > 0 && j1 < length(M) + 1 "j1 out of range"
    @assert j2 > 0 && j2 < length(M) + 1 "j2 out of range"
    @assert j1 < j2 "j1 must be less than j2"
    if j1 == 1 && j2 == length(M)
        return M
    end
    if j1 == 1
        return partial_trace_right(M, j2)
    end
    if j2 == length(M)
        return partial_trace_left(M, j1)
    else
        return partial_trace_left(partial_trace_right(M, j2), j1)
    end
end

export partial_trace