function convert_to(::Type{MPS}, v::Vector, d::Int, L::Int; cutoff=1e-12, kwargs...)
    @assert length(v) == d^L "Vector length must be d^N"
    sites = siteinds(d, L)
    T = reshape(v, (repeat([d], L)...))
    return MPS(ITensor(T, sites...), sites; cutoff=cutoff, kwargs...)
end

function convert_to(::Type{MPO}, A::Matrix, d::Int, L::Int; cutoff=1e-12, kwargs...)
    @assert size(A, 1) == d^L && size(A, 2) == d^L "Matrix dimensions must be d^N x d^N"
    sites = siteinds(d, L)
    T = reshape(A, (repeat([d], 2L)...))
    return MPO(ITensor(T, (sites..., prime.(sites)...)), sites; cutoff=cutoff, kwargs...)
end

export convert_to