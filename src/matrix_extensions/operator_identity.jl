function operator_identity(::Type{Hermitian}, n::nType) where {nType <: Integer}
    return Hermitian(Matrix{Float64}(I, n, n))
end

function operator_identity(::Type{Matrix}, n::nType) where {nType <: Integer}
    return Matrix{Float64}(I, n, n)
end

function operator_identity(::Type{SparseMatrixCSC}, n::nType) where {nType <: Integer}
    return sparse(I, n, n)
end