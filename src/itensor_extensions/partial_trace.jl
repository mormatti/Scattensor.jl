"""
    partial_trace(mpo, j1, j2)

    Computes the partial trace of an MPO over the sites between `j1` and `j2`, inclusive.
    This is obtained by contracting two by two the physical indices of the MPO to the 
    left and right of the specified positions (in the trace domain).
    
# Example
    julia> mpo = random_mpo(siteinds(3,5))
    MPO
    [1] ((dim=3|id=226|"Site,n=1")', (dim=3|id=226|"Site,n=1"), (dim=1|id=361|"Link,l=1"))
    [2] ((dim=3|id=816|"Site,n=2")', (dim=3|id=816|"Site,n=2"), (dim=1|id=505|"Link,l=2"), (dim=1|id=361|"Link,l=1"))
    [3] ((dim=3|id=229|"Site,n=3")', (dim=3|id=229|"Site,n=3"), (dim=1|id=69|"Link,l=3"), (dim=1|id=505|"Link,l=2"))
    [4] ((dim=3|id=100|"Site,n=4")', (dim=3|id=100|"Site,n=4"), (dim=1|id=862|"Link,l=4"), (dim=1|id=69|"Link,l=3"))
    [5] ((dim=3|id=357|"Site,n=5")', (dim=3|id=357|"Site,n=5"), (dim=1|id=862|"Link,l=4"))
    julia> partial_trace(mpo, 2, 4)
    MPO
    [1] ((dim=3|id=816|"Site,n=2")', (dim=3|id=816|"Site,n=2"), (dim=1|id=505|"Link,l=2"))
    [2] ((dim=3|id=229|"Site,n=3")', (dim=3|id=229|"Site,n=3"), (dim=1|id=69|"Link,l=3"), (dim=1|id=505|"Link,l=2"))
    [3] ((dim=3|id=100|"Site,n=4")', (dim=3|id=100|"Site,n=4"), (dim=1|id=69|"Link,l=3"))
"""
function partial_trace(mpo::MPO, j1::Integer, j2::Integer)
    if !(j1 > 0 && j1 < length(mpo) + 1)
        error("Position j1 must be a valid position in the MPO, got j1 = $j1")
    end
    if !(j2 > 0 && j2 < length(mpo) + 1)
        error("Position j2 must be a valid position in the MPO, got j2 = $j2")
    end
    if !(j1 <= j2)
        error("Position j1 must be less than position j2, got j1 = $j1 and j2 = $j2")
    end
    if j1 == 1 && j2 == length(mpo)
        return mpo
    end
    if j1 == 1
        return _partial_trace_right(mpo, j2)
    end
    if j2 == length(mpo)
        return _partial_trace_left(mpo, j1)
    else
        return _partial_trace_left(_partial_trace_right(mpo, j2), j1)
    end
end

function _partial_trace_left(mpo::MPO, j1::Integer)
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