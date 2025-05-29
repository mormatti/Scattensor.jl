"""
        operator_identity(MPO, d, L) -> MPO

    Returns the identity operator as an MPO (Matrix Product Operator) with bond dimension 1.
    """
function operator_identity(::Type{MPO}, d::Integer, L::Integer)
    os = OpSum()
    os += 1, "I", 1
    return truncate(MPO(os, siteinds(d, L)), cutoff = 1e-5)
end

export operator_identity