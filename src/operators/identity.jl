function operator_identity(::Type{Hermitian}, n::nType) where {nType <: Integer}
    return Hermitian(Matrix{Float64}(I, n, n))
end

function operator_identity(::Type{Matrix}, n::nType) where {nType <: Integer}
    return Matrix{Float64}(I, n, n)
end

function operator_identity(::Type{SparseMatrixCSC}, n::nType) where {nType <: Integer}
    return sparse(I, n, n)
end

function operator_identity(::Type{MPO}, d::dType, L::LType) where {dType <: Integer, LType <: Integer}
    os = OpSum()
    os += 1, "I", 1
    return truncate(MPO(os, siteinds(d, L)); cutoff = 10^(-5))
end

export operator_identity