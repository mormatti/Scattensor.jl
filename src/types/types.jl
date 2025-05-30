# Ci sei già passato. Questo è un lavoro lungo. Quando hai tempo.
abstract type HilbertSpace end

mutable struct State{HilbertSpaceType, DataType} where HilbertSpaceType <: HilbertSpace
    hilbspace::HilbertSpaceType
    data::DataType
    time::Real
end

function State{UniformChain, MPS}(wf::MPS; pbc::Bool = false)
    hilbspace = UniformChain(localdim(hs), length(hs); pbc = pbc)
    return State{UniformChain,T}(hs, wf, t)
end