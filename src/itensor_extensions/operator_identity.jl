"""
    operator_identity(MPO, d, L) -> MPO

Returns the identity operator as an MPO (Matrix Product Operator) with bond dimension 1.
The MPO is generated without link indices.

# Example
    julia> operator_identity(MPO, 3, 5)
    MPO
    [1] ((dim=3|id=830|"Site,n=1")', (dim=3|id=830|"Site,n=1"))
    [2] ((dim=3|id=79|"Site,n=2")', (dim=3|id=79|"Site,n=2"))
    [3] ((dim=3|id=902|"Site,n=3")', (dim=3|id=902|"Site,n=3"))
    [4] ((dim=3|id=839|"Site,n=4")', (dim=3|id=839|"Site,n=4"))
    [5] ((dim=3|id=302|"Site,n=5")', (dim=3|id=302|"Site,n=5"))
"""
function operator_identity(::Type{MPO}, d::Integer, L::Integer)
    s = siteinds(d, L)
    W = Vector{ITensor}(undef, L)
    for j in 1:L
        s1 = s[j]
        id_op = op("I", s1)
        W[j] = id_op
    end
    return MPO(W)
end

export operator_identity

