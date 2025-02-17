abstract type QuantumState{T} end

function wavefunction(ψ::QuantumState)
    error("Implement wavefunction for ", typeof(ψ))
end