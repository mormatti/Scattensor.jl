function dim(hs::HilbertSpace)
    error("Function dim not implemented for Hilbert space type $(typeof(hs)). Consider to implement dim(hs::$(typeof(hs))).")
end

function dim(hs::UniformChain)
    dim = (Float64(hs.localdim))^(hs.length)
    if dim == Inf
        @warn "Total Hilbert space dimension computed as infinite, consider using logdim(hs) instead."
    end
    return dim
end