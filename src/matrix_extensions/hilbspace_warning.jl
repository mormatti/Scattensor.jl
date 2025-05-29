function _hilbspace_dimension_warning(::Type{<:Matrix}, N::Real)
    if N < 0
        @warn "Invalid Hilbert space dimension, it cannot be negative. Received: N = $N"
    elseif N == Inf
        @warn "Total Hilbert space dimension of dim(H) insanely big (Infinite for floating point). Definetly not possible using dense matrices, consider using Tensor Networks."
    elseif N > 100000000
        @warn "Total Hilbert space dimension of dim(H) = $N, very large for sparse matrices, and even for sparse matrices! Consider using alternatives like Tensor Networks."
    elseif N > 20000
        @warn "Total Hilbert space dimension of dim(H) = $N, very large for dense matrices. Consider using sparse matrices."
    end
end

function _hilbspace_dimension_warning(::Type{<:SparseMatrixCSC}, N::Real)
    if N <= 0
        @warn "Invalid Hilbert space dimension, it cannot be negative. Received: N = $N"
    elseif N == Inf
        @warn "Total Hilbert space dimension of dim(H) insanely big (Infinite for floating point). Definetly not possible using sparse matrices, consider using Tensor Networks."
    elseif N > 100000000
        @warn "Total Hilbert space dimension of dim(H) = $N, very large for sparse matrices. Consider using alternatives like Tensor Networks."
    end
end

function _hilbspace_dimension_warning(::Type{<:Matrix}, d::Integer, L::Integer)
    if d < 0
        error("Invalid local dimension d, it cannot be negative. Received: d = $d")
    end
    if L < 0
        error("Invalid system size L, it cannot be negative. Received: L = $L")
    end
    if d == 0
        log10N = 0 # Not true but fine for the following warning condition
    else
        log10N = L * log10(d)
    end
    if log10N > 8
        @warn "Total Hilbert space dimension of dim(H) ∼ 10^$(Int(floor(log10N))), very large even for sparse matrices!. Consider using alternatives like Tensor Networks."
    elseif log10N > 4 + log10(2)
        @warn "Total Hilbert space dimension of dim(H) = $(d^L), very large for dense matrices. Consider using sparse matrices."
    end
end

function _hilbspace_dimension_warning(::Type{<:SparseMatrixCSC}, d::Integer, L::Integer)
    if d < 0
        error("Invalid local dimension d, it cannot be negative. Received: d = $d")
    end
    if L < 0
        error("Invalid system size L, it cannot be negative. Received: L = $L")
    end
    if d == 0
        log10N = 0 # Not true but fine for the following warning condition
    else
        log10N = L * log10(d)
    end
    if log10N > 8
        @warn "Total Hilbert space dimension of dim(H) ∼ 10^$(Int(floor(log10N))), very large for sparse matrices. Consider using alternatives like Tensor Networks."
    end
end