abstract type QuantumOperator{T} end

function representation(qo::QuantumOperator)
    error("Implement representation for ", typeof(qo))
end