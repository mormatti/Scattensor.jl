abstract type QuantumSystem end

mutable struct UniformChain <: QuantumSystem
    dims::Vector
    pbc::Bool
end