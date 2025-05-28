function operator_identity(::Type{MPO}, d::dType, L::LType) where {dType <: Integer, LType <: Integer}
    os = OpSum()
    os += 1, "I", 1
    return truncate(MPO(os, siteinds(d, L)); cutoff = 10^(-5))
end

export operator_identity