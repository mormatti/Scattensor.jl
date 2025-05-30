struct UniformChain <: HilbertSpace
    localdim::Integer
    length::Integer
    pbc::Bool
    spacing::Real
end

function UniformChain(d::Integer, L::Integer; pbc::Bool = true, a::Real = 1.0)
    return UniformChain(d, L, pbc, a)
end