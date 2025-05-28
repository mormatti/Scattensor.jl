function mps_from_vector(v::Vector, d::Int, L::Int; kwargs...)
    @assert length(v) == d^L "Vector length must be d^N"
    sites = siteinds(d, L)
    T = reshape(v, (repeat([d], L)...))
    return MPS(ITensor(T, sites...), sites; kwargs...)
end

export mps_from_vector